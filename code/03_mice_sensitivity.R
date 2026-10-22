# ============================================================
# MICE Random Forest
# 10 imputations x 10 iterations
# Parallel imputations + multithreaded ranger
# ============================================================

required_pkgs <- c(
  "here",
  "mice",
  "ranger",
  "future",
  "furrr",       # required internally by mice::futuremice()
  "parallelly"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first with:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(here)
  library(mice)
  library(ranger)
  library(future)
  library(parallelly)
})

message(">>> Running STANDALONE 03_mice_sensitivity: loads pre-MICE RDS files from disk <<<")

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

M     <- 10L
MAXIT <- 10L
NTREE <- 100L

BASE_SEED <- 123L

# Detect CPUs actually available to this R process/job
n_cores <- as.integer(parallelly::availableCores())

# Balance:
# - parallel MICE chains
# - threads available to each ranger forest
#
# At least ~2 threads/worker when enough CPUs are available.
n_workers <- min(
  M,
  max(1L, floor(n_cores / 2L))
)

threads_per_worker <- max(
  1L,
  floor(n_cores / n_workers)
)

message("Available CPU cores: ", n_cores)
message("Parallel MICE workers: ", n_workers)
message("Ranger threads per worker: ", threads_per_worker)
message(
  "Maximum concurrent threads: ~",
  n_workers * threads_per_worker
)


# ------------------------------------------------------------
# Input: load checkpoints created by 01_preprocessing.R
# ------------------------------------------------------------

data_dir <- here("analysis", "data")
mice_dir <- here("analysis", "results", "mice")

expr_path <- file.path(data_dir, "proteomic_preMICE_expr.rds")
meta_path <- file.path(data_dir, "proteomic_preMICE_metadata.rds")

message("Project root: ", here::here())
message("Loading expression checkpoint: ", expr_path)
message("Loading metadata checkpoint: ", meta_path)

if (!file.exists(expr_path)) {
  stop(
    "Missing file: ", expr_path,
    "\nRun 01_preprocessing.R first."
  )
}

if (!file.exists(meta_path)) {
  stop(
    "Missing file: ", meta_path,
    "\nRun 01_preprocessing.R first."
  )
}

expr_filtered <- readRDS(expr_path)
meta_preMICE  <- readRDS(meta_path)

expr_filtered <- as.matrix(expr_filtered)

stopifnot(
  !is.null(rownames(expr_filtered)),
  !is.null(colnames(expr_filtered)),
  "Sample_ID" %in% names(meta_preMICE),
  "sex" %in% names(meta_preMICE)
)

# Align metadata explicitly to the expression-matrix sample order.
meta_idx <- match(colnames(expr_filtered), meta_preMICE$Sample_ID)

if (anyNA(meta_idx)) {
  stop(
    "Some expression samples are absent from proteomic_preMICE_metadata.rds: ",
    paste(colnames(expr_filtered)[is.na(meta_idx)], collapse = ", ")
  )
}

meta_preMICE <- meta_preMICE[meta_idx, , drop = FALSE]

stopifnot(
  identical(
    as.character(meta_preMICE$Sample_ID),
    colnames(expr_filtered)
  )
)

# Reconstruct the same sex variable used by 01_preprocessing.R:
# sex == 1 -> male; sex == 2 -> female.
sex_raw <- meta_preMICE$sex

if (is.factor(sex_raw)) {
  sex_raw <- as.character(sex_raw)
}

sex_chr <- tolower(as.character(sex_raw))

if (all(stats::na.omit(sex_chr) %in% c("1", "2"))) {
  sex_vec <- factor(
    as.integer(sex_chr),
    levels = c(1, 2),
    labels = c("male", "female")
  )
} else if (all(stats::na.omit(sex_chr) %in% c("male", "female"))) {
  sex_vec <- factor(
    sex_chr,
    levels = c("male", "female")
  )
} else {
  stop(
    "Unexpected coding of meta_preMICE$sex. ",
    "Expected 1/2 or male/female."
  )
}

if (anyNA(sex_vec)) {
  stop("Missing or unrecognized values detected in meta_preMICE$sex.")
}

message(
  "Loaded pre-MICE data: ",
  nrow(expr_filtered), " proteins x ",
  ncol(expr_filtered), " samples"
)

message(
  "Sex counts: ",
  paste(
    names(table(sex_vec)),
    as.integer(table(sex_vec)),
    sep = "=",
    collapse = ", "
  )
)

message(
  "Missing values before MICE: ",
  sum(is.na(expr_filtered))
)

dat <- as.data.frame(
  t(expr_filtered),
  check.names = FALSE
)

stopifnot(
  nrow(dat) == length(sex_vec),
  identical(rownames(dat), colnames(expr_filtered)),
  !anyNA(sex_vec)
)

idx_by_sex <- split(
  seq_len(nrow(dat)),
  sex_vec,
  drop = TRUE
)

stopifnot(
  length(idx_by_sex) >= 1L,
  all(lengths(idx_by_sex) > 0L)
)


# ------------------------------------------------------------
# MICE separately within each sex
# ------------------------------------------------------------

imp_by_sex <- vector(
  "list",
  length(idx_by_sex)
)

names(imp_by_sex) <- names(idx_by_sex)


for (s in seq_along(idx_by_sex)) {
  
  sex_name <- names(idx_by_sex)[s]
  idx      <- idx_by_sex[[s]]
  
  message(
    "\n========================================\n",
    "MICE-RF: ", sex_name, "\n",
    "N = ", length(idx), "\n",
    "========================================"
  )
  
  dat_sex <- dat[
    idx,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Imputation method
  # Only proteins containing missing values are imputed
  # ----------------------------------------------------------
  
  meth <- rep(
    "",
    ncol(dat_sex)
  )
  
  names(meth) <- colnames(dat_sex)
  
  incomplete <- colSums(
    is.na(dat_sex)
  ) > 0
  
  meth[incomplete] <- "rf"
  
  message(
    "Proteins with missing values: ",
    sum(incomplete)
  )
  
  if (!any(incomplete)) {
    stop(
      "No missing values found in sex stratum '",
      sex_name,
      "'. futuremice() requires incomplete data."
    )
  }
  
  # A protein with no observed value in this sex stratum cannot be
  # imputed from a within-sex random-forest model.
  all_missing <- colSums(!is.na(dat_sex)) == 0L
  if (any(all_missing)) {
    stop(
      "Cannot perform within-sex MICE-RF in stratum '",
      sex_name,
      "': the following proteins are completely missing: ",
      paste(colnames(dat_sex)[all_missing], collapse = ", ")
    )
  }
  
  
  # ----------------------------------------------------------
  # Parallel MICE
  #
  # futuremice distributes the 10 imputations over workers.
  # Each RF uses ranger with its own subset of CPU threads.
  # ----------------------------------------------------------
  
  imp_by_sex[[sex_name]] <- mice::futuremice(
    data          = dat_sex,
    
    m             = M,
    maxit         = MAXIT,
    method        = meth,
    
    # Parallel MICE
    n.core        = n_workers,
    future.plan   = "multisession",
    parallelseed  = BASE_SEED + s - 1L,
    
    # Required by workers
    packages      = "ranger",
    
    # Random Forest
    rfPackage     = "ranger",
    ntree         = NTREE,
    
    # Ranger parallelism
    num.threads   = threads_per_worker,
    
    # OOB error is not needed for imputation
    oob.error     = FALSE
  )
}

# Return to the ordinary sequential future plan after parallel work.
future::plan(future::sequential)


# ------------------------------------------------------------
# Recombine men and women
# One complete dataset for each of the 10 imputations
# ------------------------------------------------------------

completed_by_m <- lapply(
  seq_len(M),
  function(m) {
    
    pieces <- lapply(
      names(idx_by_sex),
      function(sex_name) {
        
        idx <- idx_by_sex[[sex_name]]
        
        tmp <- mice::complete(
          imp_by_sex[[sex_name]],
          action = m
        )
        
        # Restore original sample IDs
        rownames(tmp) <- rownames(dat)[idx]
        
        tmp
      }
    )
    
    
    combined <- do.call(
      rbind,
      pieces
    )
    
    
    # Restore original sample and protein order
    combined <- combined[
      rownames(dat),
      colnames(dat),
      drop = FALSE
    ]
    
    combined
  }
)


# ------------------------------------------------------------
# Convert back to proteins x samples
# ------------------------------------------------------------

expr_imputed_list <- lapply(
  completed_by_m,
  function(x) {
    
    mat <- t(
      as.matrix(x)
    )
    
    rownames(mat) <- rownames(expr_filtered)
    colnames(mat) <- colnames(expr_filtered)
    
    mat
  }
)


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

na_before <- sum(
  is.na(expr_filtered)
)

na_after <- vapply(
  expr_imputed_list,
  function(x) sum(is.na(x)),
  numeric(1)
)

message("\nMissing values before MICE: ", na_before)

message(
  "Missing values after each imputation: ",
  paste(na_after, collapse = ", ")
)

stopifnot(
  length(expr_imputed_list) == M,
  all(na_after == 0)
)


# Check dimensions and names
for (i in seq_len(M)) {
  
  stopifnot(
    identical(
      dim(expr_imputed_list[[i]]),
      dim(expr_filtered)
    ),
    identical(
      rownames(expr_imputed_list[[i]]),
      rownames(expr_filtered)
    ),
    identical(
      colnames(expr_imputed_list[[i]]),
      colnames(expr_filtered)
    )
  )
}


# ------------------------------------------------------------
# Save raw MICE objects
# ------------------------------------------------------------

dir.create(
  mice_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  imp_by_sex,
  file.path(
    mice_dir,
    "mice_objects_by_sex.rds"
  )
)


saveRDS(
  meta_preMICE,
  file.path(
    mice_dir,
    "proteomic_preMICE_metadata_used.rds"
  )
)


# ------------------------------------------------------------
# Save pre-normalization imputed matrices
# ------------------------------------------------------------

for (i in seq_len(M)) {
  
  saveRDS(
    expr_imputed_list[[i]],
    file.path(
      mice_dir,
      sprintf(
        "proteomic_expr_MICE_raw_%02d.rds",
        i
      )
    )
  )
}


message(
  "\nMICE-RF completed successfully: ",
  M,
  " imputations, ",
  MAXIT,
  " iterations, ",
  NTREE,
  " trees."
)

# ============================================================
# Main-analysis robustness across MICE imputations
# ============================================================

analysis_pkgs <- c(
  "dplyr",
  "tibble",
  "ordinalNet",
  "Ckmeans.1d.dp",
  "seriation",
  "mclust",
  "limma",
  "EDASeq",
  "PhosR"
)

missing_pkgs <- analysis_pkgs[
  !vapply(analysis_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ordinalNet)
  library(Ckmeans.1d.dp)
  library(seriation)
  library(mclust)
  library(limma)
  library(EDASeq)
  library(PhosR)
})

data_dir <- here("analysis", "data")
mice_dir <- here("analysis", "results", "mice")
sens_dir <- file.path(mice_dir, "sensitivity")

dir.create(
  sens_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

M_sens <- 10L

raw_paths <- file.path(
  mice_dir,
  sprintf("proteomic_expr_MICE_raw_%02d.rds", seq_len(M_sens))
)

if (!all(file.exists(raw_paths))) {
  stop("One or more MICE expression matrices are missing.")
}

expr_mi_raw <- lapply(
  raw_paths,
  function(x) as.matrix(readRDS(x))
)

meta_pre <- readRDS(
  file.path(data_dir, "proteomic_preMICE_metadata.rds")
)

meta_raw <- readRDS(
  file.path(data_dir, "proteomic_nako_metadata_aligned.rds")
)


# ------------------------------------------------------------
# Post-MICE preprocessing
# ------------------------------------------------------------

preprocess_mice_matrix <- function(expr_imp, meta_df) {
  
  idx <- match(colnames(expr_imp), meta_df$Sample_ID)
  
  if (anyNA(idx)) {
    stop("Metadata alignment failed.")
  }
  
  meta_i <- meta_df[idx, , drop = FALSE]
  
  stopifnot(
    identical(
      as.character(meta_i$Sample_ID),
      colnames(expr_imp)
    )
  )
  
  sex_i <- factor(
    meta_i$sex,
    levels = c(1, 2),
    labels = c("male", "female")
  )
  
  age_i <- as.numeric(meta_i$age)
  
  ageclass_i <- factor(
    meta_i$ageclass,
    levels = 1:5,
    labels = c(
      "20-29",
      "30-39",
      "40-49",
      "50-59",
      "60-70"
    )
  )
  
  uq_raw <- EDASeq::betweenLaneNormalization(
    expr_imp,
    which = "upper"
  )
  
  uq_log2 <- log2(uq_raw + 1)
  
  pzp_name <- grep(
    "^(P20742|PZP(_HUMAN)?)$",
    rownames(uq_log2),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(pzp_name) != 1L) {
    stop("PZP sentinel could not be uniquely identified.")
  }
  
  pvals_by_sex <- lapply(
    levels(sex_i),
    function(sx) {
      
      ii <- which(sex_i == sx)
      
      Y <- uq_log2[, ii, drop = FALSE]
      age_s <- age_i[ii]
      
      fit <- limma::eBayes(
        limma::lmFit(
          Y,
          model.matrix(~ age_s)
        )
      )
      
      data.frame(
        protein = rownames(Y),
        pval = fit$p.value[, "age_s"]
      )
    }
  )
  
  tab <- Reduce(
    function(a, b) {
      merge(
        a,
        b,
        by = "protein",
        suffixes = c("_m", "_f")
      )
    },
    pvals_by_sex
  )
  
  tab <- tab[
    !(tab$protein %in% pzp_name),
    ,
    drop = FALSE
  ]
  
  tab$score <- pmin(
    tab$pval_m,
    tab$pval_f,
    na.rm = TRUE
  )
  
  tab <- tab[
    order(-tab$score),
    ,
    drop = FALSE
  ]
  
  if (nrow(tab) < 100L) {
    stop("Fewer than 100 negative-control candidates.")
  }
  
  controls <- head(
    tab$protein,
    100L
  )
  
  ctl_idx <- which(
    rownames(uq_log2) %in% controls
  )
  
  stopifnot(length(ctl_idx) == 100L)
  
  ruv_list <- lapply(
    levels(sex_i),
    function(sx) {
      
      ii <- which(sex_i == sx)
      
      Y <- uq_log2[, ii, drop = FALSE]
      
      M <- model.matrix(
        ~ 0 + factor(ageclass_i[ii])
      )
      
      PhosR::RUVphospho(
        Y,
        M = M,
        k = 8L,
        ctl = ctl_idx
      )
    }
  )
  
  out <- do.call(
    base::cbind,
    ruv_list
  )
  
  out <- out[
    ,
    colnames(uq_log2),
    drop = FALSE
  ]
  
  stopifnot(
    identical(dim(out), dim(expr_imp)),
    identical(rownames(out), rownames(expr_imp)),
    identical(colnames(out), colnames(expr_imp)),
    !anyNA(out)
  )
  
  out
}


expr_mi <- lapply(
  seq_len(M_sens),
  function(i) {
    
    message("Post-processing MICE dataset ", i)
    
    x <- preprocess_mice_matrix(
      expr_mi_raw[[i]],
      meta_pre
    )
    
    saveRDS(
      x,
      file.path(
        mice_dir,
        sprintf(
          "proteomic_cleaned_expr_MICE_%02d.rds",
          i
        )
      )
    )
    
    x
  }
)


# ------------------------------------------------------------
# Phenotype and primary marker sets
# ------------------------------------------------------------

sample_order <- colnames(expr_mi[[1]])

pheno_min <- meta_raw %>%
  transmute(
    sample_id = as.character(Sample_ID),
    sex = factor(
      sex,
      levels = c(1, 2),
      labels = c("male", "female")
    ),
    ageclass = factor(
      ageclass,
      levels = 1:5,
      labels = c(
        "20-29",
        "30-39",
        "40-49",
        "50-59",
        "60-70"
      )
    ),
    bmi = as.numeric(bmi),
    smok_status = factor(smok_status),
    alcohol_cons = factor(alcohol_cons),
    edu = factor(edu),
    age_years = as.numeric(age),
    age10 = as.numeric(age) / 10
  ) %>%
  filter(sample_id %in% sample_order) %>%
  arrange(match(sample_id, sample_order))

stopifnot(
  identical(
    pheno_min$sample_id,
    sample_order
  )
)

stable_ref_m <- readRDS(
  here(
    "analysis",
    "results",
    "protein_sets",
    "stable_proteins_male_ids.rds"
  )
)

stable_ref_f <- readRDS(
  here(
    "analysis",
    "results",
    "protein_sets",
    "stable_proteins_female_ids.rds"
  )
)


# ------------------------------------------------------------
# Load the exact elastic-net function from 02_analysis.R
# ------------------------------------------------------------

load_function_from_script <- function(path, object_name) {
  
  exprs <- parse(
    path,
    keep.source = FALSE
  )
  
  for (e in exprs) {
    
    is_assignment <-
      is.call(e) &&
      identical(e[[1]], as.name("<-")) &&
      is.symbol(e[[2]]) &&
      identical(
        as.character(e[[2]]),
        object_name
      )
    
    if (is_assignment) {
      eval(e, envir = .GlobalEnv)
      return(invisible(TRUE))
    }
  }
  
  stop(
    "Function not found in 02_analysis.R: ",
    object_name
  )
}

load_function_from_script(
  here("analysis", "code", "02_analysis.R"),
  "nested_enet_clust_by_sex"
)


# ------------------------------------------------------------
# Limma
# ------------------------------------------------------------

run_limma_reference <- function(
    expr_mat,
    sex_label,
    proteins
) {
  
  samples <- pheno_min$sample_id[
    pheno_min$sex == sex_label
  ]
  
  ph <- pheno_min %>%
    filter(sex == sex_label) %>%
    arrange(match(sample_id, samples))
  
  proteins <- intersect(
    proteins,
    rownames(expr_mat)
  )
  
  X <- expr_mat[
    proteins,
    samples,
    drop = FALSE
  ]
  
  keep <- is.finite(ph$age10)
  
  X <- X[
    ,
    keep,
    drop = FALSE
  ]
  
  age10 <- ph$age10[keep]
  
  design <- model.matrix(
    ~ age10
  )
  
  fit <- limma::eBayes(
    limma::lmFit(
      X,
      design
    )
  )
  
  beta <- fit$coefficients[, "age10"]
  
  se <- fit$stdev.unscaled[, "age10"] *
    sqrt(fit$s2.post)
  
  pval <- fit$p.value[, "age10"]
  
  tibble(
    protein = rownames(X),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(pval),
    fdr = p.adjust(
      pval,
      method = "BH"
    )
  )
}


# ------------------------------------------------------------
# Main analysis in each imputed dataset
# ------------------------------------------------------------

results_mi <- vector(
  "list",
  M_sens
)

for (i in seq_len(M_sens)) {
  
  message(
    "\n==============================\n",
    "Main analysis: MICE dataset ",
    i,
    "\n=============================="
  )
  
  assign(
    "expr_clean",
    expr_mi[[i]],
    envir = .GlobalEnv
  )
  
  assign(
    "pheno_min",
    pheno_min,
    envir = .GlobalEnv
  )
  
  res_m <- nested_enet_clust_by_sex(
    sex_label = "male",
    metric = "loglik",
    rule = "mean",
    nFolds_outer = 3,
    nFolds_inner = 5,
    k_clusters = 5,
    add_covariates = FALSE,
    add_interactions = FALSE,
    seed = 123
  )
  
  res_f <- nested_enet_clust_by_sex(
    sex_label = "female",
    metric = "loglik",
    rule = "mean",
    nFolds_outer = 3,
    nFolds_inner = 5,
    k_clusters = 5,
    add_covariates = FALSE,
    add_interactions = FALSE,
    seed = 123
  )
  
  stable_m <- res_m$stability %>%
    filter(selection_freq >= 2 / 3)
  
  stable_f <- res_f$stability %>%
    filter(selection_freq >= 2 / 3)
  
  limma_m <- run_limma_reference(
    expr_mi[[i]],
    "male",
    stable_ref_m
  )
  
  limma_f <- run_limma_reference(
    expr_mi[[i]],
    "female",
    stable_ref_f
  )
  
  results_mi[[i]] <- list(
    male = res_m,
    female = res_f,
    stable_m = stable_m,
    stable_f = stable_f,
    limma_m = limma_m,
    limma_f = limma_f
  )
}

saveRDS(
  results_mi,
  file.path(
    sens_dir,
    "mice_main_analysis_results.rds"
  )
)


# ------------------------------------------------------------
# Marker re-selection
# ------------------------------------------------------------

selection_summary <- function(
    results,
    sex_label,
    reference
) {
  
  sets <- lapply(
    results,
    function(x) {
      if (sex_label == "male") {
        x$stable_m$protein
      } else {
        x$stable_f$protein
      }
    }
  )
  
  all_proteins <- sort(
    unique(
      c(
        reference,
        unlist(sets)
      )
    )
  )
  
  tibble(
    protein = all_proteins,
    primary_marker = all_proteins %in% reference,
    n_MI_selected = vapply(
      all_proteins,
      function(p) {
        sum(
          vapply(
            sets,
            function(x) p %in% x,
            logical(1)
          )
        )
      },
      integer(1)
    )
  ) %>%
    mutate(
      MI_selection_frequency =
        n_MI_selected / length(sets)
    ) %>%
    arrange(
      desc(primary_marker),
      desc(MI_selection_frequency),
      protein
    )
}

selection_m <- selection_summary(
  results_mi,
  "male",
  stable_ref_m
)

selection_f <- selection_summary(
  results_mi,
  "female",
  stable_ref_f
)

write.csv(
  selection_m,
  file.path(
    sens_dir,
    "mice_marker_selection_male.csv"
  ),
  row.names = FALSE
)

write.csv(
  selection_f,
  file.path(
    sens_dir,
    "mice_marker_selection_female.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Rubin pooling and direction stability
# ------------------------------------------------------------

rubin_pool <- function(beta, se) {
  
  ok <- is.finite(beta) & is.finite(se)
  
  beta <- beta[ok]
  se <- se[ok]
  
  m <- length(beta)
  
  if (m == 0L) {
    return(
      c(
        beta = NA_real_,
        se = NA_real_,
        df = NA_real_,
        p = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }
  
  qbar <- mean(beta)
  ubar <- mean(se^2)
  
  B <- if (m > 1L) {
    stats::var(beta)
  } else {
    0
  }
  
  Tvar <- ubar +
    (1 + 1 / m) * B
  
  pooled_se <- sqrt(Tvar)
  
  if (B <= .Machine$double.eps) {
    df <- Inf
  } else {
    df <- (m - 1) *
      (
        1 +
          ubar /
          ((1 + 1 / m) * B)
      )^2
  }
  
  stat <- qbar / pooled_se
  
  pval <- if (is.finite(df)) {
    2 * stats::pt(
      -abs(stat),
      df = df
    )
  } else {
    2 * stats::pnorm(
      -abs(stat)
    )
  }
  
  crit <- if (is.finite(df)) {
    stats::qt(
      0.975,
      df = df
    )
  } else {
    stats::qnorm(0.975)
  }
  
  c(
    beta = qbar,
    se = pooled_se,
    df = df,
    p = pval,
    ci_low = qbar - crit * pooled_se,
    ci_high = qbar + crit * pooled_se
  )
}


primary_expr <- as.matrix(
  readRDS(
    file.path(
      data_dir,
      "proteomic_cleaned_expr.rds"
    )
  )
)

primary_limma_m <- run_limma_reference(
  primary_expr,
  "male",
  stable_ref_m
)

primary_limma_f <- run_limma_reference(
  primary_expr,
  "female",
  stable_ref_f
)


summarise_limma_robustness <- function(
    results,
    sex_label,
    reference,
    primary_limma
) {
  
  limma_list <- lapply(
    results,
    function(x) {
      if (sex_label == "male") {
        x$limma_m
      } else {
        x$limma_f
      }
    }
  )
  
  stable_sets <- lapply(
    results,
    function(x) {
      if (sex_label == "male") {
        x$stable_m$protein
      } else {
        x$stable_f$protein
      }
    }
  )
  
  rows <- lapply(
    reference,
    function(p) {
      
      tmp <- dplyr::bind_rows(
        lapply(
          seq_along(limma_list),
          function(i) {
            
            x <- limma_list[[i]]
            
            x <- x[
              x$protein == p,
              ,
              drop = FALSE
            ]
            
            if (nrow(x) == 0L) {
              return(NULL)
            }
            
            x$MI <- i
            x
          }
        )
      )
      
      primary_beta <- primary_limma$beta[
        match(
          p,
          primary_limma$protein
        )
      ]
      
      pooled <- rubin_pool(
        tmp$beta,
        tmp$se
      )
      
      n_selected <- sum(
        vapply(
          stable_sets,
          function(x) p %in% x,
          logical(1)
        )
      )
      
      same_direction <- mean(
        sign(tmp$beta) ==
          sign(primary_beta),
        na.rm = TRUE
      )
      
      tibble(
        protein = p,
        primary_beta = primary_beta,
        n_MI_selected = n_selected,
        MI_selection_frequency =
          n_selected / length(stable_sets),
        beta_mean = mean(
          tmp$beta,
          na.rm = TRUE
        ),
        beta_sd = sd(
          tmp$beta,
          na.rm = TRUE
        ),
        beta_min = min(
          tmp$beta,
          na.rm = TRUE
        ),
        beta_max = max(
          tmp$beta,
          na.rm = TRUE
        ),
        same_direction_frequency =
          same_direction,
        FDR05_frequency =
          mean(
            tmp$fdr < 0.05,
            na.rm = TRUE
          ),
        rubin_beta = pooled["beta"],
        rubin_se = pooled["se"],
        rubin_df = pooled["df"],
        rubin_p = pooled["p"],
        rubin_ci_low = pooled["ci_low"],
        rubin_ci_high = pooled["ci_high"]
      )
    }
  )
  
  out <- bind_rows(rows)
  
  out$rubin_fdr <- p.adjust(
    out$rubin_p,
    method = "BH"
  )
  
  out %>%
    arrange(
      desc(MI_selection_frequency),
      desc(same_direction_frequency)
    )
}

limma_robust_m <- summarise_limma_robustness(
  results_mi,
  "male",
  stable_ref_m,
  primary_limma_m
)

limma_robust_f <- summarise_limma_robustness(
  results_mi,
  "female",
  stable_ref_f,
  primary_limma_f
)

write.csv(
  limma_robust_m,
  file.path(
    sens_dir,
    "mice_limma_robustness_male.csv"
  ),
  row.names = FALSE
)

write.csv(
  limma_robust_f,
  file.path(
    sens_dir,
    "mice_limma_robustness_female.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Prediction and clustering
# ------------------------------------------------------------

confusion_metrics <- function(tab) {
  
  tab <- as.matrix(tab)
  
  total <- sum(tab)
  
  if (total == 0L) {
    return(
      tibble(
        exact_accuracy = NA_real_,
        within_one_accuracy = NA_real_
      )
    )
  }
  
  rr <- row(tab)
  cc <- col(tab)
  
  tibble(
    exact_accuracy =
      sum(diag(tab)) / total,
    within_one_accuracy =
      sum(tab[abs(rr - cc) <= 1]) / total
  )
}

performance_mi <- bind_rows(
  lapply(
    seq_along(results_mi),
    function(i) {
      
      bind_rows(
        confusion_metrics(
          results_mi[[i]]$male$confusion_overall
        ) %>%
          mutate(
            MI = i,
            sex = "male"
          ),
        
        confusion_metrics(
          results_mi[[i]]$female$confusion_overall
        ) %>%
          mutate(
            MI = i,
            sex = "female"
          )
      )
    }
  )
) %>%
  select(
    MI,
    sex,
    everything()
  )

performance_summary <- performance_mi %>%
  group_by(sex) %>%
  summarise(
    exact_mean = mean(
      exact_accuracy,
      na.rm = TRUE
    ),
    exact_sd = sd(
      exact_accuracy,
      na.rm = TRUE
    ),
    exact_min = min(
      exact_accuracy,
      na.rm = TRUE
    ),
    exact_max = max(
      exact_accuracy,
      na.rm = TRUE
    ),
    within_one_mean = mean(
      within_one_accuracy,
      na.rm = TRUE
    ),
    within_one_sd = sd(
      within_one_accuracy,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


get_clustering <- function(x, sex_label, mi) {
  
  cm <- if (sex_label == "male") {
    x$male$clustering_metrics_weighted
  } else {
    x$female$clustering_metrics_weighted
  }
  
  if (is.null(cm) || nrow(cm) == 0L) {
    return(NULL)
  }
  
  data.frame(
    MI = mi,
    sex = sex_label,
    cm,
    check.names = FALSE
  )
}

clustering_mi <- bind_rows(
  lapply(
    seq_along(results_mi),
    function(i) {
      bind_rows(
        get_clustering(
          results_mi[[i]],
          "male",
          i
        ),
        get_clustering(
          results_mi[[i]],
          "female",
          i
        )
      )
    }
  )
)

clustering_summary <- clustering_mi %>%
  group_by(sex) %>%
  summarise(
    across(
      where(is.numeric) &
        !any_of("MI"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

write.csv(
  performance_mi,
  file.path(
    sens_dir,
    "mice_prediction_by_imputation.csv"
  ),
  row.names = FALSE
)

write.csv(
  performance_summary,
  file.path(
    sens_dir,
    "mice_prediction_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  clustering_mi,
  file.path(
    sens_dir,
    "mice_clustering_by_imputation.csv"
  ),
  row.names = FALSE
)

write.csv(
  clustering_summary,
  file.path(
    sens_dir,
    "mice_clustering_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# AgeProt
# ------------------------------------------------------------

ageprot_r2 <- function(
    expr_mat,
    sex_label,
    stable_ids
) {
  
  samples <- pheno_min$sample_id[
    pheno_min$sex == sex_label
  ]
  
  ph <- pheno_min %>%
    filter(sex == sex_label) %>%
    arrange(match(sample_id, samples))
  
  proteins <- intersect(
    stable_ids,
    rownames(expr_mat)
  )
  
  X <- expr_mat[
    proteins,
    samples,
    drop = FALSE
  ]
  
  keep <- is.finite(ph$age10)
  
  X <- X[
    ,
    keep,
    drop = FALSE
  ]
  
  ph <- ph[
    keep,
    ,
    drop = FALSE
  ]
  
  X_z <- t(
    scale(
      t(X),
      center = TRUE,
      scale = TRUE
    )
  )
  
  X_z[!is.finite(X_z)] <- 0
  
  design <- model.matrix(
    ~ age10,
    data = ph
  )
  
  fit <- limma::eBayes(
    limma::lmFit(
      X_z,
      design
    )
  )
  
  beta <- fit$coefficients[, "age10"]
  beta[!is.finite(beta)] <- 0
  
  score <- as.numeric(
    t(X_z) %*% beta
  )
  
  fit_score <- lm(
    score ~ ph$age10
  )
  
  tibble(
    n = length(score),
    n_proteins = length(proteins),
    R2 = summary(fit_score)$r.squared
  )
}

primary_ageprot_m <- ageprot_r2(
  primary_expr,
  "male",
  stable_ref_m
)$R2

primary_ageprot_f <- ageprot_r2(
  primary_expr,
  "female",
  stable_ref_f
)$R2

ageprot_mi <- bind_rows(
  lapply(
    seq_len(M_sens),
    function(i) {
      
      bind_rows(
        ageprot_r2(
          expr_mi[[i]],
          "male",
          stable_ref_m
        ) %>%
          mutate(
            MI = i,
            sex = "male",
            primary_R2 = primary_ageprot_m
          ),
        
        ageprot_r2(
          expr_mi[[i]],
          "female",
          stable_ref_f
        ) %>%
          mutate(
            MI = i,
            sex = "female",
            primary_R2 = primary_ageprot_f
          )
      )
    }
  )
) %>%
  mutate(
    delta_R2 = R2 - primary_R2
  ) %>%
  select(
    MI,
    sex,
    everything()
  )

ageprot_summary <- ageprot_mi %>%
  dplyr::group_by(sex) %>%
  dplyr::summarise(
    primary_R2 = primary_R2[[1L]],
    MICE_R2_mean = mean(
      R2,
      na.rm = TRUE
    ),
    MICE_R2_sd = stats::sd(
      R2,
      na.rm = TRUE
    ),
    MICE_R2_min = min(
      R2,
      na.rm = TRUE
    ),
    MICE_R2_max = max(
      R2,
      na.rm = TRUE
    ),
    mean_delta_R2 = mean(
      delta_R2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write.csv(
  ageprot_mi,
  file.path(
    sens_dir,
    "mice_ageprot_by_imputation.csv"
  ),
  row.names = FALSE
)

write.csv(
  ageprot_summary,
  file.path(
    sens_dir,
    "mice_ageprot_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\n========================================\n",
  "MICE SENSITIVITY ANALYSIS\n",
  "========================================\n"
)

cat("\nPrimary male markers:\n")
print(
  selection_m %>%
    dplyr::filter(primary_marker)
)

cat("\nPrimary female markers:\n")
print(
  selection_f %>%
    dplyr::filter(primary_marker)
)

cat("\nMale limma robustness:\n")
print(
  limma_robust_m %>%
    dplyr::select(
      protein,
      MI_selection_frequency,
      same_direction_frequency,
      FDR05_frequency,
      primary_beta,
      rubin_beta,
      rubin_ci_low,
      rubin_ci_high,
      rubin_p,
      rubin_fdr
    )
)

cat("\nFemale limma robustness:\n")
print(
  limma_robust_f %>%
    dplyr::select(
      protein,
      MI_selection_frequency,
      same_direction_frequency,
      FDR05_frequency,
      primary_beta,
      rubin_beta,
      rubin_ci_low,
      rubin_ci_high,
      rubin_p,
      rubin_fdr
    )
)

cat("\nPrediction robustness:\n")
print(performance_summary)

cat("\nClustering robustness:\n")
print(clustering_summary)

cat("\nAgeProt robustness:\n")
print(ageprot_summary)

capture.output(
  {
    cat(
      "MICE sensitivity analysis\n",
      "=========================\n\n"
    )
    
    cat("Male primary marker re-selection\n")
    print(
      selection_m %>%
        dplyr::filter(primary_marker)
    )
    
    cat("\nFemale primary marker re-selection\n")
    print(
      selection_f %>%
        dplyr::filter(primary_marker)
    )
    
    cat("\nMale limma robustness\n")
    print(limma_robust_m)
    
    cat("\nFemale limma robustness\n")
    print(limma_robust_f)
    
    cat("\nPrediction robustness\n")
    print(performance_summary)
    
    cat("\nClustering robustness\n")
    print(clustering_summary)
    
    cat("\nAgeProt robustness\n")
    print(ageprot_summary)
  },
  file = file.path(
    sens_dir,
    "mice_sensitivity_summary.txt"
  )
)

message(
  "\nSensitivity analysis completed. Results saved in: ",
  sens_dir
)
