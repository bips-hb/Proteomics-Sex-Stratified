# 02_analysis.R
# Main analyses and manuscript outputs.
# Run after 01_preprocessing.R.


suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(tibble)
  library(purrr)
  library(ordinalNet)
  library(missForest)
  library(Ckmeans.1d.dp)
  library(seriation)
  library(mclust)
  library(readxl)
  library(limma)
  library(scales)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggpattern)
  library(ggtree)
  library(ggnewscale)
  library(grid)
  library(gridExtra)
  library(pheatmap)
  library(gt)
  library(bootnet)
  library(qgraph)
  library(igraph)
})

here::i_am("analysis/code/02_analysis.R")

set.seed(123)

data_dir <- here("analysis", "data")
results_dir <- here("analysis", "results")
fig_main_dir <- file.path(results_dir, "figures_main")
fig_supp_dir <- file.path(results_dir, "figures_supplementary")
tables_dir <- file.path(results_dir, "tables")
protein_dir <- file.path(results_dir, "protein_sets")
network_dir <- file.path(results_dir, "bootnet_aging_driven")

invisible(lapply(
  c(fig_main_dir, fig_supp_dir, tables_dir, protein_dir, network_dir),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# Typography is shared by all main and supplementary figures.
PLOT_FONT <- "Arial"
for (.pkg in c("ragg", "systemfonts")) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop("Install the graphics package '", .pkg, "' before running this script.")
  }
}
if (!any(tolower(systemfonts::system_fonts()$family) == tolower(PLOT_FONT))) {
  stop("Arial is not installed. Install it, restart R, and run again.")
}
if (!capabilities("cairo")) {
  stop("Cairo support is required to export PDF figures with Arial.")
}
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts("Arial" = grDevices::windowsFont(PLOT_FONT))
}
ggplot2::theme_update(text = ggplot2::element_text(family = PLOT_FONT))
ggplot2::update_geom_defaults("text", list(family = PLOT_FONT))
ggplot2::update_geom_defaults("label", list(family = PLOT_FONT))
ggplot2::update_geom_defaults(ggrepel::GeomTextRepel, list(family = PLOT_FONT))
ggplot2::update_geom_defaults(ggrepel::GeomLabelRepel, list(family = PLOT_FONT))


.save_figure_device <- function(draw_fun, destination, width_px, height_px,
                                dpi, format, pointsize = 12) {
  tmp <- tempfile(pattern = "nako_figure_", fileext = paste0(".", format))
  on.exit(unlink(tmp), add = TRUE)
  if (format == "png") {
    ragg::agg_png(tmp, width = width_px, height = height_px, units = "px",
                  res = dpi, pointsize = pointsize, background = "white")
  } else {
    grDevices::cairo_pdf(tmp, width = width_px / dpi, height = height_px / dpi,
                         family = PLOT_FONT, pointsize = pointsize, bg = "white")
  }
  device_id <- grDevices::dev.cur()
  tryCatch({
    graphics::par(family = PLOT_FONT)
    draw_fun()
  }, finally = {
    if (device_id %in% grDevices::dev.list()) grDevices::dev.off(device_id)
  })
  if (!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size <= 0) {
    stop("The graphics device did not complete: ", destination)
  }
  if (!file.copy(tmp, destination, overwrite = TRUE) ||
      !file.exists(destination) || file.info(destination)$size != file.info(tmp)$size) {
    stop("Could not save the completed figure to: ", destination)
  }
  invisible(destination)
}

save_base_plot <- function(draw_fun, stem, width_px, height_px, dpi = 96) {
  for (format in c("png", "pdf")) {
    .save_figure_device(draw_fun, paste0(stem, ".", format),
                        width_px, height_px, dpi, format)
  }
  invisible(c(png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf")))
}

save_gg <- function(p, stem, width_px, height_px, dpi = 96) {
  save_base_plot(function() print(p), stem, width_px, height_px, dpi)
  invisible(p)
}

save_main_base <- function(draw_fun, stem, height_mm, width_mm = 180,
                           dpi = 600, pointsize = 8) {
  width_px <- width_mm / 25.4 * dpi
  height_px <- height_mm / 25.4 * dpi
  for (format in c("png", "pdf")) {
    .save_figure_device(draw_fun, paste0(stem, ".", format),
                        width_px, height_px, dpi, format, pointsize)
  }
  invisible(c(png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf")))
}

save_main_gg <- function(p, stem, height_mm, width_mm = 180, dpi = 600) {
  save_main_base(function() print(p), stem, height_mm, width_mm, dpi)
  invisible(p)
}

save_gt <- function(tbl, stem, vwidth = 1952, vheight = 900) {
  # Always save formats that do not require an interactive browser dependency.
  suppressMessages(
    gt::gtsave(tbl, paste0(stem, ".html"))
  )
  
  try(
    suppressMessages(
      gt::gtsave(tbl, paste0(stem, ".rtf"))
    ),
    silent = TRUE
  )
  
  # PNG export requires webshot2. Never let gt open an interactive
  # installation menu while the analysis script is running.
  if (requireNamespace("webshot2", quietly = TRUE)) {
    suppressMessages(
      gt::gtsave(
        tbl,
        paste0(stem, ".png"),
        vwidth = vwidth,
        vheight = vheight,
        expand = 0
      )
    )
  }
  
  invisible(tbl)
}

# Publication labels only: keep original protein IDs internally
PROTEIN_LABEL_MAP <- c(
  ALS = "IGFALS", ANGT = "AGT", ATPA = "ATP5F1A", ATS13 = "ADAMTS13",
  CERU = "CP", CO6A1 = "COL6A1", GPIX = "GP9", HV307 = "IGHV3-7",
  HV315 = "IGHV3-15", IBP2 = "IGFBP2", IGJ = "JCHAIN", ITB3 = "ITGB3",
  LG3BP = "LGALS3BP", LV211 = "IGLV2-11", LV469 = "IGLV4-69",
  LYSC = "LYZ", PAFA = "PLA2G7", PROP = "CFP", TENX = "TNXB", TRFL = "LTF"
)

pub_protein <- function(x) {
  y <- sub("_HUMAN$", "", as.character(x), ignore.case = TRUE)
  z <- unname(PROTEIN_LABEL_MAP[toupper(y)])
  z[is.na(z)] <- y[is.na(z)]
  z
}

# 1. Load processed data and phenotype
expr_path <- here("analysis", "data", "proteomic_cleaned_expr.rds")
meta_path <- here("analysis", "data", "proteomic_nako_metadata_aligned.rds")

expr_clean <- readRDS(expr_path)   # matrix: proteins x samples (post-filtering, imputation, UQ, RUV)
meta_raw   <- readRDS(meta_path)   # aligned metadata (same samples as in preprocessing)

# Quick checks
stopifnot(is.matrix(expr_clean) || is.data.frame(expr_clean))
expr_clean <- as.matrix(expr_clean)
cat("Expression matrix dims:", paste(dim(expr_clean), collapse = " x "), "\n")
cat("Metadata rows:", nrow(meta_raw), " | columns:", ncol(meta_raw), "\n")

pheno_min <- meta_raw %>%
  transmute(
    sample_id   = Sample_ID,
    sex         = factor(sex, levels = c(1, 2), labels = c("male", "female")),
    ageclass    = factor(ageclass, levels = 1:5,
                         labels = c("20-29", "30-39", "40-49", "50-59", "60-70")),
    bmi         = as.numeric(bmi),
    smok_status = factor(smok_status),
    alcohol_cons = factor(alcohol_cons),
    edu         = factor(edu),
    age_years   = as.numeric(age)
  )

# Keep only samples present in the expression matrix (defensive)
pheno_min <- pheno_min %>%
  filter(sample_id %in% colnames(expr_clean)) %>%
  arrange(match(sample_id, colnames(expr_clean)))

# Final alignment check
stopifnot(identical(pheno_min$sample_id, colnames(expr_clean)))

cat("\nPhenotype summary:\n")
print(summary(pheno_min %>% select(-sample_id)))
cat("\nSex counts:\n"); print(table(pheno_min$sex, useNA = "ifany"))
cat("Ageclass counts:\n"); print(table(pheno_min$ageclass, useNA = "ifany"))

# Impute 'edu' with missForest within each sex stratum

imp_frame <- pheno_min %>%
  dplyr::select(sample_id, sex, ageclass, bmi, smok_status, alcohol_cons, edu)

by_sex <- split(imp_frame, imp_frame$sex)

imp_list <- lapply(by_sex, function(df_sex) {
  sids <- df_sex$sample_id
  x <- df_sex %>% dplyr::select(-sample_id, -sex)
  
  set.seed(123)
  fit <- missForest(x, ntree = 500, verbose = FALSE)
  
  data.frame(
    sample_id = sids,
    sex       = df_sex$sex,
    edu_imp   = fit$ximp$edu,
    row.names = NULL
  )
})

imp_by_sex <- dplyr::bind_rows(imp_list)

imp_by_sex <- imp_by_sex %>%
  dplyr::arrange(match(sample_id, pheno_min$sample_id))
stopifnot(identical(imp_by_sex$sample_id, pheno_min$sample_id))

orig_levels <- levels(pheno_min$edu)
pheno_min$edu <- factor(imp_by_sex$edu_imp, levels = union(orig_levels, unique(imp_by_sex$edu_imp)))

cat("NAs in edu after imputation:", sum(is.na(pheno_min$edu)), "\n")
print(table(pheno_min$sex, pheno_min$edu, useNA = "ifany"))

# 2. Ordinal elastic net + held-out unsupervised validation
nested_enet_clust_by_sex <- function(
    sex_label,
    alpha = 0.5,
    metric = c("loglik","devPct","misclass","brier"),
    rule   = c("mean","1se"),
    nFolds_outer = 5,
    nFolds_inner = 5,
    k_clusters   = 5,
    add_covariates   = TRUE,
    add_interactions = FALSE,
    maxiterOut = 50,
    seed = 123
) {
  ## ------------------------------------------------------------
  ## prerequisites: expr_clean (proteins x samples), pheno_min aligned
  ## ------------------------------------------------------------
  stopifnot(exists("expr_clean"), exists("pheno_min"))
  metric <- match.arg(metric)
  rule   <- match.arg(rule)
  set.seed(seed)
  
  ## ---------- helpers ----------
  # folds balanced with age
  make_stratified_folds <- function(y, K) {
    y <- droplevels(y)
    lv <- levels(y)
    fold <- integer(length(y))
    for (lev in lv) {
      idx <- which(y == lev)
      idx <- sample(idx)
      fold[idx] <- rep(1:K, length.out = length(idx))
    }
    fold
  }
  
  pick_lambda <- function(tunefit, metric, rule, nFolds) {
    Z <- switch(metric,
                loglik   = tunefit$loglik,
                devPct   = tunefit$devPct,
                misclass = -tunefit$misclass,
                brier    = -tunefit$brier
    )
    m  <- rowMeans(Z, na.rm = TRUE)
    se <- apply(Z, 1, function(v) stats::sd(v, na.rm = TRUE)) / sqrt(nFolds)
    best <- which.max(m)
    if (rule == "1se") {
      thr <- m[best] - se[best]
      ok  <- which(m >= thr)
      best <- ok[which.max(tunefit$lambdaVals[ok])]
    }
    best
  }
  
  # clustering metrics
  purity_score <- function(tab) sum(apply(tab, 1, max)) / sum(tab)
  
  .adjacency_index <- function(cluster, age_idx) {
    df <- data.frame(cluster = as.integer(cluster), age = age_idx)
    by(df$age, df$cluster, function(v) {
      n <- length(v); if (n < 2) return(NA_real_)
      tab <- table(v); lev <- as.integer(names(tab))
      total <- choose(n, 2); good <- 0
      for (i in seq_along(lev)) {
        a <- lev[i]; na <- tab[i]
        good <- good + choose(na, 2)
        if ((a + 1) %in% lev) good <- good + na * tab[as.character(a + 1)]
      }
      good / total
    }) |> unlist() |> mean(na.rm = TRUE)
  }
  .span_score <- function(cluster, age_idx) {
    df <- data.frame(cluster = as.integer(cluster), age = age_idx)
    tabs <- split(df$age, df$cluster)
    spans <- vapply(tabs, function(v) max(v) - min(v), numeric(1))
    w <- vapply(tabs, length, integer(1))
    sum(spans * w) / sum(w)
  }
  .cluster_rank_by_mean_age <- function(cluster, age_idx) {
    mu <- tapply(age_idx, cluster, mean)
    ord <- order(mu)
    as.integer(factor(cluster, levels = names(mu)[ord]))
  }
  .ordinal_spearman <- function(cluster, age_idx) {
    cr <- .cluster_rank_by_mean_age(cluster, age_idx)
    ct <- suppressWarnings(cor.test(age_idx, cr, method = "spearman", exact = FALSE))
    list(rho = unname(ct$estimate), p.value = unname(ct$p.value))
  }
  .kendall_tau_global <- function(cluster, age_idx) {
    cr <- .cluster_rank_by_mean_age(cluster, age_idx)
    suppressWarnings(cor(age_idx, cr, method = "kendall"))
  }
  .between_cluster_concordance <- function(cluster, age_idx) {
    cr <- .cluster_rank_by_mean_age(cluster, age_idx)
    n <- length(cr); if (n < 2) return(NA_real_)
    C <- utils::combn(n, 2L)
    dc <- cr[C[2,]] - cr[C[1,]]; da <- age_idx[C[2,]] - age_idx[C[1,]]
    sel <- (dc != 0) & (da != 0); if (!any(sel)) return(NA_real_)
    mean(sign(dc[sel]) == sign(da[sel]))
  }
  
  run_ordered_clustering_seriation <- function(X_prot, age_fac = NULL,
                                               pc_select   = "elbow",
                                               kmax_elbow  = 30L,
                                               order_method= "OLO",
                                               k_mode      = "fixed",
                                               k_fixed     = 5L,
                                               seed = 123) {
    stopifnot(is.matrix(X_prot))
    set.seed(seed)
    
    # samples x features
    Xs <- t(X_prot)
    keep_var <- apply(Xs, 2, function(v) stats::var(v, na.rm = TRUE) > 0)
    Xs <- Xs[, keep_var, drop = FALSE]
    
    for (j in seq_len(ncol(Xs))) {
      if (anyNA(Xs[, j])) {
        m <- mean(Xs[, j], na.rm = TRUE)
        Xs[is.na(Xs[, j]), j] <- m
      }
    }
    Xs <- scale(Xs, center = TRUE, scale = TRUE)
    
    pca <- prcomp(Xs, center = FALSE, scale. = FALSE)
    eig <- pca$sdev^2
    if (pc_select == "elbow") {
      v <- eig; v <- v[is.finite(v) & v > 0]
      kmax_elbow <- max(2L, min(kmax_elbow, length(v)))
      x <- seq_len(kmax_elbow); y <- v[seq_len(kmax_elbow)]
      num <- abs((y[kmax_elbow]-y[1])*x - (kmax_elbow-1)*y + kmax_elbow*y[1] - y[kmax_elbow])
      den <- sqrt((y[kmax_elbow]-y[1])^2 + (kmax_elbow-1)^2)
      q <- which.max(num/den)
    } else {
      cum <- cumsum(eig)/sum(eig)
      q <- which(cum >= 0.90)[1]
    }
    q <- max(1L, as.integer(q))
    S <- pca$x[, seq_len(q), drop = FALSE]
    
    D <- dist(S)
    ser_ctrl <- if (order_method == "OLO") list(linkage = "ward.D2") else list()
    ord <- seriation::get_order(seriation::seriate(D, method = order_method, control = ser_ctrl))
    
    Dm <- as.matrix(D)
    path_increments <- c(0, Dm[cbind(ord[-length(ord)], ord[-1])])
    z <- cumsum(path_increments)
    z <- as.numeric(scale(z))
    
    if (k_mode == "fixed") {
      k_opt <- max(2L, as.integer(k_fixed))
    } else {
      stop("only k_mode='fixed' implemented here")
    }
    
    cl_raw <- Ckmeans.1d.dp::Ckmeans.1d.dp(z, k_opt)$cluster
    ordlabs <- order(tapply(z, cl_raw, mean))
    cl <- as.integer(factor(cl_raw, levels = unique(cl_raw)[ordlabs]))
    cl <- factor(cl, levels = sort(unique(cl)))
    
    conf <- NULL; ari <- pur <- rho <- rho_p <- tau <- adj <- span <- bcc <- NA_real_
    if (!is.null(age_fac)) {
      age_fac <- droplevels(age_fac)
      age_idx <- as.integer(age_fac)
      conf <- table(Cluster = cl, AgeGroup = age_fac)
      ari  <- mclust::adjustedRandIndex(cl, age_fac)
      pur  <- purity_score(conf)
      tmp  <- .ordinal_spearman(cl, age_idx)
      rho  <- tmp$rho; rho_p <- tmp$p.value
      tau  <- .kendall_tau_global(cl, age_idx)
      adj  <- .between_cluster_concordance(cl, age_idx)
      span <- .span_score(cl, age_idx)
      bcc  <- .between_cluster_concordance(cl, age_idx)
    }
    
    list(
      cluster = cl,
      PCs_used = q,
      confusion = conf,
      ARI = ari, Purity = pur,
      OrdinalRho = rho, Rho_pvalue = rho_p,
      KendallTau = tau, AdjIndex = adj, Span = span, BCC = bcc
    )
  }
  
  ## ---------- 1) sex subset and build X,y ----------
  sel_samples <- pheno_min$sample_id[pheno_min$sex == sex_label]
  sel_samples <- intersect(sel_samples, colnames(expr_clean))
  Xprot_all <- t(expr_clean[, sel_samples, drop = FALSE])
  
  meta_sel <- pheno_min |>
    dplyr::filter(sample_id %in% sel_samples) |>
    dplyr::arrange(match(sample_id, sel_samples))
  
  y_all <- ordered(meta_sel$ageclass, levels = levels(pheno_min$ageclass))
  
  Xfull_all <- Xprot_all
  if (add_covariates) {
    cov_df <- meta_sel[, c("bmi","smok_status","alcohol_cons","edu")]
    cov_df$bmi <- as.numeric(cov_df$bmi)
    if (anyNA(cov_df$bmi)) cov_df$bmi[is.na(cov_df$bmi)] <- median(cov_df$bmi, na.rm = TRUE)
    to_missing <- function(f) factor(ifelse(is.na(f), "Missing", as.character(f)))
    cov_df$smok_status  <- to_missing(cov_df$smok_status)
    cov_df$alcohol_cons <- to_missing(cov_df$alcohol_cons)
    cov_df$edu          <- to_missing(cov_df$edu)
    Xcov <- model.matrix(~ scale(bmi) + smok_status + alcohol_cons + edu - 1, data = cov_df)
    Xfull_all <- cbind(Xfull_all, Xcov)
  }
  if (add_interactions) {
    Xfull_all <- stats::model.matrix(~ .^2 - 1, data = as.data.frame(Xfull_all))
  }
  
  keep <- complete.cases(y_all)
  Xfull_all <- Xfull_all[keep, , drop = FALSE]
  y_all     <- y_all[keep]
  samples_all <- sel_samples[keep]
  
  ## ---------- 2) outer folds ----------
  folds <- make_stratified_folds(y_all, nFolds_outer)
  
  fold_results <- vector("list", nFolds_outer)
  agg_conf <- NULL
  
  for (k in seq_len(nFolds_outer)) {
    tr_idx <- which(folds != k)
    te_idx <- which(folds == k)
    
    y_tr <- y_all[tr_idx]
    y_te <- y_all[te_idx]
    
    X_tr <- Xfull_all[tr_idx, , drop = FALSE]
    X_te_samples <- samples_all[te_idx]
    
    if (length(unique(y_tr)) < 2L || length(te_idx) == 0L) {
      warning(sprintf("[%s] fold %d skipped (too few samples/classes)", sex_label, k))
      fold_results[[k]] <- list(
        fold = k,
        selected_proteins = character(0),
        coef_table = NULL,
        clustering = NULL,
        test_samples = X_te_samples
      )
      next
    }
    
    # --- inner CV on training ---
    tune_tr <- ordinalNet::ordinalNetTune(
      x = X_tr, y = y_tr,
      family = "cumulative",
      parallelTerms = TRUE,
      nonparallelTerms = FALSE,
      alpha = alpha,
      nFolds = nFolds_inner,
      standardize = TRUE,
      maxiterOut = maxiterOut,
      printProgress = FALSE
    )
    best_idx <- pick_lambda(tune_tr, metric, rule, nFolds_inner)
    
    # --- coefficients at chosen λ (trained on TRAIN only!) ---
    B <- coef(tune_tr$fit, whichLambda = best_idx, matrix = TRUE)
    B_df <- B |>
      as.data.frame() |>
      tibble::rownames_to_column("predictor") |>
      dplyr::filter(predictor != "(Intercept)")
    
    # summarize
    cmat <- as.matrix(B_df[, -1, drop = FALSE])
    B_df$beta_mean   <- rowMeans(cmat)
    B_df$beta_maxabs <- apply(abs(cmat), 1, max)
    
    # selected vars in TRAIN
    sel_mask <- apply(cmat, 1, function(v) any(abs(v) > 0))
    coef_table_tr <- B_df[sel_mask, , drop = FALSE] |>
      dplyr::arrange(dplyr::desc(beta_maxabs))
    
    # --- keep ONLY proteins for clustering ---
    prot_selected <- intersect(coef_table_tr$predictor, rownames(expr_clean))
    
    if (length(prot_selected) >= 2) {
      X_test_prot <- expr_clean[prot_selected, X_te_samples, drop = FALSE]
      cl_res <- run_ordered_clustering_seriation(
        X_prot = X_test_prot,
        age_fac = y_te,
        pc_select = "elbow",
        kmax_elbow = 30L,
        order_method = "OLO",
        k_mode = "fixed",
        k_fixed = k_clusters,
        seed = seed + k
      )
      # aggregate confusion
      if (!is.null(cl_res$confusion)) {
        if (is.null(agg_conf)) {
          agg_conf <- cl_res$confusion
        } else {
          all_r <- union(rownames(agg_conf), rownames(cl_res$confusion))
          all_c <- union(colnames(agg_conf), colnames(cl_res$confusion))
          tmp1 <- matrix(0, nrow = length(all_r), ncol = length(all_c),
                         dimnames = list(all_r, all_c))
          tmp2 <- tmp1
          tmp1[rownames(agg_conf), colnames(agg_conf)] <- agg_conf
          tmp2[rownames(cl_res$confusion), colnames(cl_res$confusion)] <- cl_res$confusion
          agg_conf <- tmp1 + tmp2
        }
      }
    } else {
      cl_res <- list(
        cluster = factor(rep(1L, length(y_te))),
        PCs_used = NA_integer_,
        confusion = NULL,
        ARI = NA_real_, Purity = NA_real_,
        OrdinalRho = NA_real_, Rho_pvalue = NA_real_,
        KendallTau = NA_real_, AdjIndex = NA_real_,
        Span = NA_real_, BCC = NA_real_
      )
    }
    
    fold_results[[k]] <- list(
      fold = k,
      lambda_index = best_idx,
      selected_proteins = prot_selected,
      coef_table = coef_table_tr,
      clustering = cl_res,
      test_samples = X_te_samples
    )
  }
  
  ## ---------- 3) stability selection across folds ----------
  all_prots <- sort(unique(unlist(lapply(fold_results, function(fr) fr$selected_proteins))))
  if (length(all_prots)) {
    beta_mat <- sapply(fold_results, function(fr) {
      if (is.null(fr$coef_table)) return(rep(NA_real_, length(all_prots)))
      mm <- match(all_prots, fr$coef_table$predictor)
      out <- rep(NA_real_, length(all_prots))
      out[!is.na(mm)] <- fr$coef_table$beta_mean[mm[!is.na(mm)]]
      out
    })
    if (!is.matrix(beta_mat)) beta_mat <- matrix(beta_mat, nrow = length(all_prots))
    sel_freq <- rowMeans(!is.na(beta_mat))
    mean_beta <- apply(beta_mat, 1, function(v) mean(v, na.rm = TRUE))
    median_beta <- apply(beta_mat, 1, function(v) median(v, na.rm = TRUE))
    stability_tbl <- data.frame(
      protein = all_prots,
      selection_freq = sel_freq,
      mean_beta = mean_beta,
      median_beta = median_beta,
      stringsAsFactors = FALSE
    ) |>
      dplyr::arrange(dplyr::desc(selection_freq), dplyr::desc(abs(mean_beta)))
  } else {
    stability_tbl <- data.frame(
      protein = character(0),
      selection_freq = numeric(0),
      mean_beta = numeric(0),
      median_beta = numeric(0)
    )
  }
  
  ## ---------- 4) average clustering metrics ----------
  cl_rows <- lapply(seq_along(fold_results), function(i) {
    fr <- fold_results[[i]]
    if (is.null(fr$clustering)) return(NULL)
    n_test <- length(fr$test_samples)
    data.frame(
      fold = i,
      n_test = n_test,
      ARI = fr$clustering$ARI,
      Purity = fr$clustering$Purity,
      OrdinalRho = fr$clustering$OrdinalRho,
      KendallTau = fr$clustering$KendallTau,
      AdjIndex = fr$clustering$AdjIndex,
      Span = fr$clustering$Span,
      BCC = fr$clustering$BCC
    )
  })
  cl_rows <- dplyr::bind_rows(cl_rows)
  if (nrow(cl_rows)) {
    w <- cl_rows$n_test / sum(cl_rows$n_test)
    agg_metrics <- with(cl_rows, data.frame(
      ARI = sum(ARI * w, na.rm = TRUE),
      Purity = sum(Purity * w, na.rm = TRUE),
      OrdinalRho = sum(OrdinalRho * w, na.rm = TRUE),
      KendallTau = sum(KendallTau * w, na.rm = TRUE),
      AdjIndex = sum(AdjIndex * w, na.rm = TRUE),
      Span = sum(Span * w, na.rm = TRUE),
      BCC = sum(BCC * w, na.rm = TRUE)
    ))
  } else {
    agg_metrics <- NULL
  }
  
  ## ---------- 5) return ----------
  invisible(list(
    sex = sex_label,
    outer_folds = nFolds_outer,
    clustering_metrics_weighted = agg_metrics,
    confusion_overall = agg_conf,
    stability = stability_tbl,
    folds = fold_results
  ))
}

# 3 min
res_male_nested <- nested_enet_clust_by_sex(
  sex_label      = "male",
  metric         = "loglik",
  rule           = "mean",
  nFolds_outer   = 3,
  nFolds_inner   = 5,
  k_clusters     = 5,
  add_covariates = FALSE,
  add_interactions = FALSE,
  seed = 123
)

res_male_nested$clustering_metrics_weighted
res_male_nested$confusion_overall


prot_names_path <- here("analysis", "data", "0009_BIPS_NAKO_ProteinTable.xlsx")
prot_info <- read_excel(prot_names_path, sheet = 1)
prot_info <- prot_info[-1, ] %>%
  select(PG.ProteinDescriptions, PG.ProteinNames) %>%
  dplyr::rename(description = PG.ProteinDescriptions,
                protein = PG.ProteinNames) %>%
  distinct()

stable_prots_m <- res_male_nested$stability %>%
  filter(selection_freq >= 2/3)

# The sign of the limma coefficient corresponds to the sign of the marginal association with age.

sel_samples_m <- pheno_min$sample_id[pheno_min$sex == "male"]
stable_proteins_m <- stable_prots_m$protein
X_limma_m <- expr_clean[stable_proteins_m, sel_samples_m, drop = FALSE]


age_years_m <- pheno_min$age_years[pheno_min$sex == "male"]
age10_m <- age_years_m / 10

keep_m <- is.finite(age10_m)
X_limma_m <- X_limma_m[, keep_m, drop = FALSE]
age10_m <- age10_m[keep_m]

design_m <- model.matrix(~ age10_m)
stopifnot(nrow(design_m) == ncol(X_limma_m))

fit_m <- lmFit(X_limma_m, design_m)
fit_m <- eBayes(fit_m)

beta_m <- fit_m$coefficients[, "age10_m"]
pval_m <- fit_m$p.value[, "age10_m"]
fdr_m  <- p.adjust(pval_m, method = "BH")

limma_m <- tibble::tibble(
  protein  = rownames(X_limma_m),
  beta_age = as.numeric(beta_m),     
  p_age    = as.numeric(pval_m),
  fdr_age  = as.numeric(fdr_m)
) %>%
  dplyr::left_join(
    stable_prots_m %>% dplyr::select(protein, selection_freq),
    by = "protein"
  ) %>%
  mutate(selection_freq = round(selection_freq, 2)) %>%
  dplyr::arrange(dplyr::desc(beta_age))   

stable_annot_m <- limma_m %>%
  left_join(prot_info, by = "protein")

stable_annot_m

res_female_nested <- nested_enet_clust_by_sex(
  sex_label      = "female",
  metric         = "loglik",
  rule           = "mean",
  nFolds_outer   = 3,
  nFolds_inner   = 5,
  k_clusters     = 5,
  add_covariates = FALSE,
  add_interactions = FALSE,
  seed = 123
)

res_female_nested$clustering_metrics_weighted
res_female_nested$confusion_overall

stable_prots_f <- res_female_nested$stability %>%
  filter(selection_freq >= 2/3)

# The sign of the limma coefficient corresponds to the sign of the marginal association with age.

sel_samples_f <- pheno_min$sample_id[pheno_min$sex == "female"]
stable_proteins_f <- stable_prots_f$protein
X_limma_f <- expr_clean[stable_proteins_f, sel_samples_f, drop = FALSE]


age_years_f <- pheno_min$age_years[pheno_min$sex == "female"]
age10_f <- age_years_f / 10

keep_f <- is.finite(age10_f)
X_limma_f <- X_limma_f[, keep_f, drop = FALSE]
age10_f <- age10_f[keep_f]

design_f <- model.matrix(~ age10_f)
stopifnot(nrow(design_f) == ncol(X_limma_f))

fit_f <- lmFit(X_limma_f, design_f)
fit_f <- eBayes(fit_f)

beta_f <- fit_f$coefficients[, "age10_f"]
pval_f <- fit_f$p.value[, "age10_f"]
fdr_f  <- p.adjust(pval_f, method = "BH")

limma_f <- tibble::tibble(
  protein  = rownames(X_limma_f),
  beta_age = as.numeric(beta_f),     
  p_age    = as.numeric(pval_f),
  fdr_age  = as.numeric(fdr_f)
) %>%
  dplyr::left_join(
    stable_prots_f %>% dplyr::select(protein, selection_freq),
    by = "protein"
  ) %>%
  dplyr::mutate(selection_freq = round(selection_freq, 2)) %>%
  dplyr::arrange(dplyr::desc(beta_age))   

stable_annot_f <- limma_f %>%
  dplyr::left_join(prot_info, by = "protein")

stable_annot_f

out_dir <- here("analysis", "results", "protein_sets")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

saveRDS(
  stable_annot_m,
  file = file.path(out_dir, "stable_proteins_male_annotated.rds")
)

saveRDS(
  stable_annot_f,
  file = file.path(out_dir, "stable_proteins_female_annotated.rds")
)

saveRDS(
  stable_annot_m$protein,
  file = file.path(out_dir, "stable_proteins_male_ids.rds")
)

saveRDS(
  stable_annot_f$protein,
  file = file.path(out_dir, "stable_proteins_female_ids.rds")
)

cat("File salvati in:\n", out_dir, "\n")

# 3. Table S2
library(dplyr)
library(gt)

# Helpers
# ----------------------------
get_age_coef_name <- function(fit) {
  cn <- colnames(fit$coefficients)
  cand <- cn[grepl("^age10", cn) | grepl("^age", cn)]
  cand <- setdiff(cand, "(Intercept)")
  if (length(cand) >= 1) return(cand[1])
  if (length(cn) == 2) return(setdiff(cn, "(Intercept)")[1])
  stop("Impossibile identificare il coefficiente età in fit$coefficients. Colonne: ", paste(cn, collapse = ", "))
}

limma_ci95 <- function(fit, coef_name) {
  stopifnot(coef_name %in% colnames(fit$coefficients))
  beta <- fit$coefficients[, coef_name]
  se   <- fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post)
  tcr  <- stats::qt(0.975, df = fit$df.total)
  
  tibble(
    protein = rownames(fit$coefficients),
    ci_low  = as.numeric(beta - tcr * se),
    ci_high = as.numeric(beta + tcr * se)
  )
}

fmt_ci <- function(lo, hi) sprintf("%.2f to %.2f", lo, hi)

make_sex_block <- function(limma_df, fit_obj, sex_label,
                           beta_col = "beta_age",
                           p_col = "p_age",
                           fdr_col = "fdr_age",
                           sel_col = "selection_freq",
                           stable_cut = 2/3) {
  
  stopifnot(all(c("protein", beta_col, p_col, fdr_col, sel_col) %in% names(limma_df)))
  
  coef_name <- get_age_coef_name(fit_obj)
  ci_df <- limma_ci95(fit_obj, coef_name)
  
  limma_df %>%
    inner_join(ci_df, by = "protein") %>%
    filter(.data[[sel_col]] >= stable_cut) %>%
    transmute(
      protein = pub_protein(protein),
      beta    = round(.data[[beta_col]], 2),
      ci_95   = fmt_ci(ci_low, ci_high),
      p_raw   = .data[[p_col]],
      fdr     = .data[[fdr_col]],
      sel     = if_else(abs(.data[[sel_col]] - 1) < 1e-8, "3/3", "2/3")
    ) %>%
    arrange(desc(beta)) %>%
    mutate(rank = row_number()) %>%
    select(rank, everything())
}

# ----------------------------
# 1) Build blocks (Male / Female) from your objects
# ----------------------------
male_block   <- make_sex_block(limma_m, fit_m, "Male")
female_block <- make_sex_block(limma_f, fit_f, "Female")

# Prefix columns (except rank) for spanners
male_w <- male_block %>%
  rename_with(~ paste0("Male_", .x), -rank)

female_w <- female_block %>%
  rename_with(~ paste0("Female_", .x), -rank)

tab_df <- full_join(male_w, female_w, by = "rank") %>%
  arrange(rank)

male_cols   <- c("Male_protein","Male_beta","Male_ci_95","Male_p_raw","Male_fdr","Male_sel")
female_cols <- c("Female_protein","Female_beta","Female_ci_95","Female_p_raw","Female_fdr","Female_sel")

# ----------------------------
# 2) GT table (single table, clean)
# ----------------------------
gt_tbl <- tab_df %>%
  gt() %>%
  tab_header(title = "Selected age-associated proteins, stratified by sex") %>%
  cols_label(
    rank = "",
    Male_protein = "Protein", Male_beta = "Beta", Male_ci_95 = "95% CI", Male_p_raw = "p", Male_fdr = "FDR", Male_sel = "Selection",
    Female_protein = "Protein", Female_beta = "Beta", Female_ci_95 = "95% CI", Female_p_raw = "p", Female_fdr = "FDR", Female_sel = "Selection"
  ) %>%
  tab_spanner(label = "Male",   columns = all_of(male_cols)) %>%
  tab_spanner(label = "Female", columns = all_of(female_cols)) %>%
  fmt_number(columns = c(Male_p_raw, Male_fdr, Female_p_raw, Female_fdr), decimals = 3) %>%
  sub_missing(everything(), missing_text = "") %>%
  cols_align(align = "left",   columns = c(Male_protein, Female_protein)) %>%
  cols_align(align = "center", columns = c(setdiff(male_cols, "Male_protein"),
                                           setdiff(female_cols, "Female_protein"))) %>%
  tab_options(
    table.width = pct(100),
    table.font.size = px(10),
    data_row.padding = px(2)
  ) %>%
  tab_footnote(
    footnote = "Beta, 95% CI, p-values, and FDR are from sex-stratified limma models (effect per 10 years). Selection indicates stability across outer folds (2/3 or 3/3).",
    locations = cells_title(groups = "title")
  )

gt_tbl

# Optional save:
# output is saved below
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# gt::gtsave(gt_tbl, filename = file.path(out_dir, "Table2_limma_coeffs.html"))
# gt::gtsave(gt_tbl, filename = file.path(out_dir, "Table2_limma_coeffs.png"), expand = 0)
# try(gt::gtsave(gt_tbl, filename = file.path(out_dir, "Table2_limma_coeffs.rtf")))

save_gt(
  gt_tbl,
  file.path(tables_dir, "Table2_limma_coeffs"),
  vwidth = 1952, vheight = 850
)

# 4. Figure 2: volcano + limma coefficient bars
## ============================================================
## COMBINED FIGURE:
##   Panel A = Volcano (Male/Female) + custom legend
##   Panel B = Barplot (Male/Female) with stable-only limma betas
## Requirements (assumed to exist from earlier pipeline):
##   - expr_clean (proteins x samples)
##   - pheno_min  (aligned; has sample_id, sex, age_years)
##   - res_male_nested, res_female_nested (from nested_enet_clust_by_sex)
##   - prot_info (optional; used only if label_mode_use = "description")
## Packages: dplyr, ggplot2, patchwork, limma, ggrepel, ggpattern, scales
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(limma)
  library(ggrepel)
  library(ggpattern)
  library(scales)
})

stopifnot(exists("expr_clean"), exists("pheno_min"),
          exists("res_male_nested"), exists("res_female_nested"))

## ============================================================
## 1) VOLCANO: limma on ALL proteins + fold stability annotation
## ============================================================

run_limma_all <- function(expr, pheno, sex_label) {
  stopifnot(all(c("sample_id","sex","age_years") %in% names(pheno)))
  
  keep_sex <- (pheno$sex == sex_label) & (pheno$sample_id %in% colnames(expr))
  sel_samples <- pheno$sample_id[keep_sex]
  X <- expr[, sel_samples, drop = FALSE]
  
  age10 <- as.numeric(pheno$age_years[keep_sex]) / 10
  keep_age <- is.finite(age10)
  
  X <- X[, keep_age, drop = FALSE]
  age10 <- age10[keep_age]
  
  design <- model.matrix(~ age10)
  stopifnot(nrow(design) == ncol(X))
  
  fit <- limma::lmFit(X, design)
  fit <- limma::eBayes(fit)
  
  tibble(
    protein = rownames(X),
    beta    = as.numeric(fit$coefficients[, "age10"]),
    pval    = as.numeric(fit$p.value[, "age10"]),
    fdr     = as.numeric(p.adjust(fit$p.value[, "age10"], method = "BH"))
  )
}

annotate_enet_folds_sign <- function(limma_df, stability_tbl, nFolds_outer) {
  stopifnot(all(c("protein","selection_freq") %in% names(stability_tbl)))
  stopifnot(is.numeric(nFolds_outer) && nFolds_outer >= 2)
  
  stab <- stability_tbl %>%
    select(protein, selection_freq) %>%
    mutate(
      folds_selected = as.integer(round(selection_freq * nFolds_outer)),
      folds_selected = pmax(pmin(folds_selected, nFolds_outer), 0L)
    )
  
  limma_df %>%
    left_join(stab, by = "protein") %>%
    mutate(
      folds_selected = if_else(is.na(folds_selected), 0L, folds_selected),
      fold_group = case_when(
        folds_selected >= nFolds_outer ~ paste0(nFolds_outer, "/", nFolds_outer, " folds"),
        folds_selected == (nFolds_outer - 1L) ~ paste0((nFolds_outer - 1L), "/", nFolds_outer, " folds"),
        TRUE ~ paste0("0–", (nFolds_outer - 2L), "/", nFolds_outer, " folds")
      ),
      # for colors we keep only the 2 top strata explicitly, everything else grey
      color_group = case_when(
        folds_selected == (nFolds_outer - 1L) & beta > 0 ~ "Positive (top-1)",
        folds_selected == (nFolds_outer - 1L) & beta < 0 ~ "Negative (top-1)",
        folds_selected >= nFolds_outer       & beta > 0 ~ "Positive (top)",
        folds_selected >= nFolds_outer       & beta < 0 ~ "Negative (top)",
        TRUE ~ "Other"
      ),
      color_group = factor(
        color_group,
        levels = c("Positive (top-1)","Negative (top-1)","Positive (top)","Negative (top)","Other")
      )
    )
}

plot_volcano_folds_sign <- function(df,
                                    title,
                                    p_nominal = 0.1,
                                    col_other = "#BDBDBD",
                                    col_pos_top1 = "#FF8FD1",
                                    col_neg_top1 = "#2EC5FF",
                                    col_pos_top  = "#FF0033",
                                    col_neg_top  = "#0B4DFF",
                                    pt_other  = 1.2,
                                    pt_sel    = 1.8,
                                    label_size = 2.5,
                                    base_size  = 8,
                                    axis_text  = 7.5,
                                    axis_title = 9,
                                    title_size = 10) {
  
  stopifnot(is.finite(p_nominal), p_nominal > 0, p_nominal < 1)
  y_thr <- -log10(p_nominal)
  
  col_map <- c(
    "Positive (top-1)" = col_pos_top1,
    "Negative (top-1)" = col_neg_top1,
    "Positive (top)"   = col_pos_top,
    "Negative (top)"   = col_neg_top,
    "Other"            = col_other
  )
  
  df_lab <- df %>%
    filter(color_group %in% c("Positive (top-1)","Negative (top-1)","Positive (top)","Negative (top)")) %>%
    mutate(label = pub_protein(protein))
  
  ggplot(df, aes(x = beta, y = -log10(pval))) +
    geom_point(
      data = df %>% filter(color_group == "Other"),
      aes(color = color_group),
      size = pt_other,
      alpha = 0.25
    ) +
    geom_point(
      data = df %>% filter(color_group != "Other"),
      aes(color = color_group),
      size = pt_sel,
      alpha = 1
    ) +
    ggrepel::geom_text_repel(
      data = df_lab,
      aes(label = label, color = color_group),
      size = label_size,
      min.segment.length = 0,
      family = PLOT_FONT,
      box.padding = 0.30,
      point.padding = 0.25,
      point.size = ifelse(df_lab$label == "PRDX6", pt_sel, 1),
      nudge_x = ifelse(df_lab$label == "PRDX6",
                       -0.08 * diff(range(df$beta, na.rm = TRUE)), 0),
      force = 2,
      max.iter = 100000,
      max.time = 5,
      max.overlaps = Inf,
      seed = 123,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.35, color = "black") +
    scale_color_manual(values = col_map, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
    labs(
      title = title,
      x = "limma coefficient (beta per 10 years)",
      y = expression(-log[10](p))
    ) +
    theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    theme(
      plot.title   = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      axis.text.x  = element_text(family = PLOT_FONT, size = axis_text),
      axis.text.y  = element_text(family = PLOT_FONT, size = axis_text),
      axis.title.x = element_text(family = PLOT_FONT, size = axis_title),
      axis.title.y = element_text(family = PLOT_FONT, size = axis_title),
      plot.margin  = margin(3, 5, 3, 5)
    )
}

plot_volcano_legend_1row6col <- function(col_grey_m = "#9EBCAD",
                                         col_grey_f = "#CBB69A",
                                         col_pos2   = "#FF8FD1",
                                         col_neg2   = "#2EC5FF",
                                         col_pos3   = "#FF0033",
                                         col_neg3   = "#0B4DFF",
                                         base_size = 8,
                                         text_size = 7,
                                         # layout controls
                                         x_start = 1.2,
                                         dx = 1.55) {
  
  title_text <- "Selection frequency and direction:"
  
  pts <- tibble::tibble(
    x = x_start + dx * (0:5),
    y = 1,
    col = c(col_grey_m, col_grey_f, col_pos2, col_neg2, col_pos3, col_neg3),
    alpha = c(0.6, 0.6, 1, 1, 1, 1),
    lab = c("0–1/3 folds M", "0–1/3 folds F",
            "2/3 folds (β>0)", "2/3 folds (β<0)",
            "3/3 folds (β>0)", "3/3 folds (β<0)")
  )
  
  x_end <- max(pts$x) + dx
  
  ggplot2::ggplot() +
    # titolo sopra
    ggplot2::annotate("text",
                      x = x_start, y = 1.28, label = title_text,
                      hjust = 0, vjust = 0.5, fontface = "bold",
                      size = text_size / ggplot2::.pt
    ) +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = x, y = y, colour = col, alpha = alpha),
      shape = 16, size = 1.8, show.legend = FALSE
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::geom_text(
      data = pts,
      ggplot2::aes(x = x, y = y, label = lab),
      hjust = 0, vjust = 0.5,
      nudge_x = 0.10,
      size = text_size / ggplot2::.pt
    ) +
    ggplot2::coord_cartesian(xlim = c(x_start, x_end), ylim = c(0.85, 1.40), clip = "off") +
    ggplot2::theme_void(base_size = base_size, base_family = PLOT_FONT) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 5, 0, 5))
}



limma_all_m <- run_limma_all(expr_clean, pheno_min, "male")
limma_all_f <- run_limma_all(expr_clean, pheno_min, "female")

volc_m <- annotate_enet_folds_sign(limma_all_m, res_male_nested$stability, res_male_nested$outer_folds)
volc_f <- annotate_enet_folds_sign(limma_all_f, res_female_nested$stability, res_female_nested$outer_folds)

p_volc_m <- plot_volcano_folds_sign(volc_m, "Male", p_nominal = 0.1, col_other = "#9EBCAD")
p_volc_f <- plot_volcano_folds_sign(volc_f, "Female", p_nominal = 0.1, col_other = "#CBB69A")

p_leg <- plot_volcano_legend_1row6col()

p_volcano_block <-
  ((p_volc_m + p_volc_f) + plot_layout(ncol = 2, widths = c(1, 1))) / p_leg +
  plot_layout(heights = c(1, 0.10))

## ============================================================
## 2) BARPLOT: limma on STABLE proteins only + pattern by stability
## ============================================================

run_limma_stable_only <- function(expr_mat, pheno, sex_label, stable_tbl) {
  
  sel_samples <- intersect(
    pheno$sample_id[pheno$sex == sex_label],
    colnames(expr_mat)
  )
  
  stable_prots <- intersect(stable_tbl$protein, rownames(expr_mat))
  stopifnot(length(stable_prots) > 0)
  
  X <- expr_mat[stable_prots, sel_samples, drop = FALSE]
  
  meta_sel <- pheno %>%
    filter(sex == sex_label, sample_id %in% sel_samples) %>%
    arrange(match(sample_id, sel_samples))
  
  age10 <- as.numeric(meta_sel$age_years) / 10
  
  keep <- is.finite(age10)
  X <- X[, keep, drop = FALSE]
  age10 <- age10[keep]
  
  design <- model.matrix(~ age10)
  stopifnot(nrow(design) == ncol(X))
  
  fit <- limma::lmFit(X, design)
  fit <- limma::eBayes(fit)
  
  tibble(
    protein = rownames(X),
    beta    = as.numeric(fit$coefficients[, "age10"]),
    pval    = as.numeric(fit$p.value[, "age10"]),
    fdr     = as.numeric(p.adjust(fit$p.value[, "age10"], method = "BH"))
  )
}

prep_stable_bar_beta_df <- function(stability_tbl,
                                    limma_stable,
                                    prot_info = NULL,
                                    min_freq = 2/3,
                                    top_n = 30,
                                    label_mode = c("protein", "description"),
                                    nFolds_outer = 3) {
  
  label_mode <- match.arg(label_mode)
  
  df <- stability_tbl %>%
    select(protein, selection_freq) %>%
    inner_join(limma_stable %>% select(protein, beta), by = "protein") %>%
    filter(selection_freq >= min_freq)
  
  if (!is.null(prot_info)) {
    df <- df %>% left_join(prot_info %>% select(protein, description), by = "protein")
  } else {
    df$description <- NA_character_
  }
  
  df <- df %>%
    mutate(
      protein_short = pub_protein(protein),
      label = case_when(
        label_mode == "description" & !is.na(description) & description != "" ~ description,
        TRUE ~ protein_short
      )
    )
  
  # Map selection_freq -> folds label robustly for generic nFolds
  eps <- 1e-8
  df <- df %>%
    mutate(
      folds_selected = as.integer(round(selection_freq * nFolds_outer)),
      freq_lab = case_when(
        folds_selected == (nFolds_outer - 1L) ~ paste0((nFolds_outer - 1L), "/", nFolds_outer),
        folds_selected >= nFolds_outer        ~ paste0(nFolds_outer, "/", nFolds_outer),
        TRUE ~ NA_character_
      ),
      freq_lab = factor(freq_lab, levels = c(paste0((nFolds_outer - 1L), "/", nFolds_outer),
                                             paste0(nFolds_outer, "/", nFolds_outer)))
    ) %>%
    filter(!is.na(freq_lab))
  
  df %>%
    arrange(desc(beta)) %>%
    slice_head(n = top_n) %>%
    mutate(label = factor(label, levels = rev(unique(label))))
}

plot_stable_bar_beta <- function(df,
                                 title,
                                 sex_col,
                                 base_fill = "grey85",
                                 base_size = 8,
                                 axis_text = 7,
                                 axis_title = 8,
                                 title_size = 10,
                                 x_pad = 0.10) {
  
  if (nrow(df) == 0) {
    return(ggplot() + theme_void(base_family = PLOT_FONT) + ggtitle(paste0(title, ": no proteins after filter")))
  }
  
  xmax <- max(abs(df$beta), na.rm = TRUE)
  xmax <- xmax * (1 + x_pad)
  
  levs <- levels(df$freq_lab)
  stopifnot(length(levs) == 2)
  
  pattern_vals <- setNames(c("stripe", "none"), levs)
  fill_vals    <- setNames(c(base_fill, sex_col), levs)
  lab_vals     <- setNames(paste0(levs, " folds"), levs)
  
  ggplot(df, aes(y = label, x = beta, fill = freq_lab, pattern = freq_lab)) +
    ggpattern::geom_col_pattern(
      color = "grey25",
      linewidth = 0.15,
      width = 0.8,
      pattern_fill   = sex_col,
      pattern_color  = sex_col,
      pattern_density = 0.55,
      pattern_spacing = 0.03,
      pattern_key_scale_factor = 0.7
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, alpha = 0.7) +
    scale_pattern_manual(
      values = pattern_vals,
      breaks = levs,
      labels = lab_vals,
      name   = "Selection frequency"
    ) +
    scale_fill_manual(
      values = fill_vals,
      breaks = levs,
      labels = lab_vals,
      name   = "Selection frequency"
    ) +
    scale_x_continuous(
      limits = c(-xmax, xmax),
      labels = function(x) sprintf("%.2f", x),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(title = title, x = "limma coefficient (beta per 10 years)", y = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    theme(
      plot.title   = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      axis.text.x  = element_text(family = PLOT_FONT, size = axis_text),
      axis.text.y  = element_text(family = PLOT_FONT, size = axis_text),
      axis.title.x = element_text(family = PLOT_FONT, size = axis_title),
      legend.position = "bottom",
      legend.key.size = unit(3.4, "mm"),
      legend.text = element_text(family = PLOT_FONT, size = axis_text),
      legend.title = element_text(family = PLOT_FONT, size = axis_text, face = 'bold'),
      plot.margin = margin(3, 5, 3, 4)
    )
}

# ---- stable sets (>= 2/3) ----
stable_prots_m <- res_male_nested$stability %>% filter(selection_freq >= 2/3)
stable_prots_f <- res_female_nested$stability %>% filter(selection_freq >= 2/3)

limma_stable_m <- run_limma_stable_only(expr_clean, pheno_min, "male", stable_prots_m)
limma_stable_f <- run_limma_stable_only(expr_clean, pheno_min, "female", stable_prots_f)

label_mode_use <- "protein"   # or "description"
top_n_use      <- 30

df_m_bar <- prep_stable_bar_beta_df(
  stability_tbl  = res_male_nested$stability,
  limma_stable   = limma_stable_m,
  prot_info      = if (exists("prot_info")) prot_info else NULL,
  min_freq       = 2/3,
  top_n          = top_n_use,
  label_mode     = label_mode_use,
  nFolds_outer   = res_male_nested$outer_folds
)

df_f_bar <- prep_stable_bar_beta_df(
  stability_tbl  = res_female_nested$stability,
  limma_stable   = limma_stable_f,
  prot_info      = if (exists("prot_info")) prot_info else NULL,
  min_freq       = 2/3,
  top_n          = top_n_use,
  label_mode     = label_mode_use,
  nFolds_outer   = res_female_nested$outer_folds
)

p_bar_m <- plot_stable_bar_beta(df_m_bar, "Male",   sex_col = "#009E73")
p_bar_f <- plot_stable_bar_beta(df_f_bar, "Female", sex_col = "#E69F00")

p_bar_block <- (p_bar_m + p_bar_f) + plot_layout(ncol = 2, widths = c(1, 1))

## ============================================================
## 3) FINAL COMBINATION:
##   - volcano block on top, tagged "A"
##   - barplot block below, tagged "B"
##   - NO A/B tags inside the subpanels
## ============================================================



p_final_combined <-
  (wrap_elements(full = p_volcano_block) /
     wrap_elements(full = p_bar_block)) +
  plot_layout(heights = c(1.05, 1.00)) +
  plot_annotation(tag_levels = "A") &   
  theme(
    plot.tag = element_text(family = PLOT_FONT, size = 12, face = "bold"),
    plot.tag.position = c(0.01, 0.99),  
    plot.margin = margin(3, 5, 3, 5)
  )

p_final_combined

save_main_gg(
  p_final_combined,
  file.path(fig_main_dir, "volcano_barplot"),
  height_mm = 180
)

# 5. Figure 3: confusion matrices + circular dendrograms
## ============================================================
## COMBINED FIGURE (2 BLOCKS):
##   A) Confusion matrices (Male/Female)  [TOP]
##   B) Circular tree + 2 rings (Male/Female + legends) [BOTTOM]
##
## Requirements already in your env:
##   - expr_clean, pheno_min
##   - res_male_nested, res_female_nested
##   - stable_prots_m, stable_prots_f  (>= 2/3)  [or we rebuild them below]
##
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggtree)
  library(limma)
  library(rlang)
  library(grid)
  library(ggnewscale)
})

## ============================================================
## 0) Defensive objects
## ============================================================
stopifnot(exists("expr_clean"), exists("pheno_min"))
# 6. Supplementary Figures S3-S5
stopifnot(exists("res_male_nested"), exists("res_female_nested"))

# Rebuild stable sets if not present
if (!exists("stable_prots_m")) {
  stable_prots_m <- res_male_nested$stability %>% dplyr::filter(selection_freq >= 2/3)
}
if (!exists("stable_prots_f")) {
  stable_prots_f <- res_female_nested$stability %>% dplyr::filter(selection_freq >= 2/3)
}

## ============================================================
## 1) BLOCK A (TOP): confusion matrices (Male / Female)
## ============================================================

conf_to_df <- function(M) {
  M <- as.matrix(M)
  stopifnot(!is.null(rownames(M)), !is.null(colnames(M)))
  
  df <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  names(df) <- c("AgeGroup", "Cluster", "n")
  
  df$AgeGroup <- factor(df$AgeGroup, levels = rownames(M))
  df$Cluster  <- factor(df$Cluster,  levels = colnames(M))
  df
}

plot_confmat_panel <- function(M,
                               title = NULL,
                               palette_color = "#009E73",
                               show_counts = TRUE,
                               digits_perc = 0,
                               normalize = c("none","row","col"),
                               base_size   = 8,
                               axis_text   = 7.5,
                               axis_title  = 8.5,
                               title_size  = 10) {
  
  normalize <- match.arg(normalize)
  df <- conf_to_df(M)
  
  if (normalize == "row") {
    df <- df %>% group_by(AgeGroup) %>% mutate(val = n / sum(n)) %>% ungroup()
    fill_lab <- "Row %"
    lab <- if (show_counts) paste0(df$n, "\n(", percent(df$val, accuracy = 10^-digits_perc), ")") else
      percent(df$val, accuracy = 10^-digits_perc)
  } else if (normalize == "col") {
    df <- df %>% group_by(Cluster) %>% mutate(val = n / sum(n)) %>% ungroup()
    fill_lab <- "Col %"
    lab <- if (show_counts) paste0(df$n, "\n(", percent(df$val, accuracy = 10^-digits_perc), ")") else
      percent(df$val, accuracy = 10^-digits_perc)
  } else {
    df$val <- df$n
    fill_lab <- "Count"
    lab <- if (show_counts) as.character(df$n) else ""
  }
  
  df <- df %>%
    mutate(
      diag = as.integer(AgeGroup) == as.integer(Cluster),
      lab  = lab
    )
  
  ggplot(df, aes(x = Cluster, y = AgeGroup, fill = val)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_tile(data = df %>% filter(diag),
              fill = NA, color = "black", linewidth = 0.5) +
    geom_text(aes(label = lab), size = 2.5) +
    scale_fill_gradient(
      low = "white",
      high = palette_color,
      labels = if (normalize == "none") comma else percent_format(),
      name = fill_lab
    ) +
    scale_y_discrete(
      limits = rev,
      labels = c(
        "20-29" = "20-29",
        "30-39" = "30-39",
        "40-49" = "40-49",
        "50-59" = "50-59",
        "60-70" = "60-69"
      )
    ) +
    labs(title = title, x = "Cluster (predicted)", y = "Age group (true)") +
    theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    theme(
      plot.title   = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      axis.text.x  = element_text(family = PLOT_FONT, size = axis_text),
      axis.text.y  = element_text(family = PLOT_FONT, size = axis_text),
      axis.title.x = element_text(family = PLOT_FONT, size = axis_title),
      axis.title.y = element_text(family = PLOT_FONT, size = axis_title),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.width = unit(3, "mm"),
      legend.key.height = unit(4, "mm"),
      plot.margin = margin(3, 4, 3, 4)
    )
}

# confusion matrices (AgeGroup x Cluster)
conf_m <- t(res_male_nested$confusion_overall)
conf_f <- t(res_female_nested$confusion_overall)

max_count <- max(conf_m, conf_f, na.rm = TRUE)

p_conf_m <- plot_confmat_panel(
  conf_m, title = "Male",
  palette_color = "#009E73",
  normalize = "none",
  show_counts = TRUE
) + scale_fill_gradient(low = "white", high = "#009E73",
                        limits = c(0, max_count), name = "Count")

p_conf_f <- plot_confmat_panel(
  conf_f, title = "Female",
  palette_color = "#E69F00",
  normalize = "none",
  show_counts = TRUE
) + scale_fill_gradient(low = "white", high = "#E69F00",
                        limits = c(0, max_count), name = "Count")

p_conf_block <-
  (p_conf_m + p_conf_f) +
  plot_layout(ncol = 2, widths = c(1, 1)) &
  theme(
    plot.margin   = margin(3, 5, 3, 5),
    panel.spacing = unit(4, "mm")
  )

## ============================================================
## 2) BLOCK B (BOTTOM): circular tree (Male/Female) + legends
## ============================================================

## limma on stable only (ring beta)
run_limma_stable_only <- function(expr_mat, pheno, sex_label, stable_tbl) {
  
  sel_samples <- intersect(
    pheno$sample_id[pheno$sex == sex_label],
    colnames(expr_mat)
  )
  
  stable_prots <- intersect(stable_tbl$protein, rownames(expr_mat))
  stopifnot(length(stable_prots) > 0)
  
  X <- expr_mat[stable_prots, sel_samples, drop = FALSE]
  
  meta_sel <- pheno %>%
    filter(sex == sex_label, sample_id %in% sel_samples) %>%
    arrange(match(sample_id, sel_samples))
  
  age10 <- as.numeric(meta_sel$age_years) / 10
  keep <- is.finite(age10)
  X <- X[, keep, drop = FALSE]
  age10 <- age10[keep]
  
  design <- model.matrix(~ age10)
  stopifnot(nrow(design) == ncol(X))
  
  fit <- limma::lmFit(X, design)
  fit <- limma::eBayes(fit)
  
  tibble(
    protein = rownames(X),
    beta    = as.numeric(fit$coefficients[, "age10"]),
    pval    = as.numeric(fit$p.value[, "age10"]),
    fdr     = as.numeric(p.adjust(fit$p.value[, "age10"], method = "BH"))
  )
}

# (re)fit if needed
limma_stable_m <- run_limma_stable_only(expr_clean, pheno_min, "male",   stable_prots_m)
limma_stable_f <- run_limma_stable_only(expr_clean, pheno_min, "female", stable_prots_f)

SEX_COL <- c(male = "#009E73", female = "#E69F00")

tint_with_white <- function(col, t = 0.65) {
  rgb <- grDevices::col2rgb(col) / 255
  out <- rgb * (1 - t) + 1 * t
  grDevices::rgb(out[1], out[2], out[3])
}

get_cluster_cols <- function(sex_label) {
  sex_col <- SEX_COL[[sex_label]]
  c(
    "1" = tint_with_white(sex_col, t = 0.70),
    "2" = sex_col
  )
}

.strip_human <- function(x) sub("_HUMAN$", "", x)

.build_tree_base <- function(hc, title = NULL, base_size = 8) {
  p <- ggtree(hc, layout = "circular", size = 0.25) +
    theme_void(base_size = base_size, base_family = PLOT_FONT) +
    theme(plot.margin = margin(3, 3, 3, 3)) +
    labs(title = title)
  
  tips <- p$data %>%
    dplyr::filter(.data$isTip) %>%
    dplyr::arrange(.data$y) %>%
    dplyr::mutate(
      angle = 360 * (row_number() - 0.5) / n(),
      hjust = ifelse(angle > 90 & angle < 270, 1, 0),
      angle2 = ifelse(angle > 90 & angle < 270, angle + 180, angle)
    )
  
  list(plot = p, tips = tips)
}

.add_ring <- function(p_tree, tips_df, ring_df, value_col,
                      ring_inner, ring_thick,
                      fill_scale,
                      legend_title,
                      show_legend = TRUE) {
  
  df <- tips_df %>%
    dplyr::select(label, x, y, angle, hjust, angle2) %>%
    dplyr::left_join(ring_df, by = "label") %>%
    dplyr::mutate(
      x_ring = ring_inner + 0.5 * ring_thick,
      w = ring_thick,
      h = 1.00
    )
  
  p_tree +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = df,
      aes(x = x_ring, y = y, width = w, height = h, fill = .data[[value_col]]),
      color = "grey70", linewidth = 0.15
    ) +
    fill_scale +
    guides(fill = if (show_legend) guide_legend(title = legend_title) else "none")
}

.add_tip_labels <- function(p_tree, tips_df, r_label, label_size = 3.2) {
  df_lab <- tips_df %>% dplyr::mutate(label2 = pub_protein(label))
  p_tree +
    geom_text(
      data = df_lab,
      aes(x = r_label, y = y, label = label2, angle = angle2, hjust = hjust),
      size = label_size
    )
}

build_circular_panel_2rings <- function(
    sex_label = c("male","female"),
    hc_method = "average",
    label_size = 3.2
) {
  sex_label <- match.arg(sex_label)
  
  if (sex_label == "male") {
    limma_st  <- limma_stable_m
    stability <- stable_prots_m
  } else {
    limma_st  <- limma_stable_f
    stability <- stable_prots_f
  }
  
  prot <- intersect(stability$protein, limma_st$protein)
  prot <- prot[!is.na(prot)]
  stopifnot(length(prot) >= 2)
  
  sel_samples <- pheno_min$sample_id[pheno_min$sex == sex_label]
  sel_samples <- intersect(sel_samples, colnames(expr_clean))
  X <- expr_clean[prot, sel_samples, drop = FALSE]
  
  C <- suppressWarnings(cor(t(X), use = "pairwise.complete.obs", method = "pearson"))
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  
  d  <- as.dist(1 - C)
  hc <- hclust(d, method = hc_method)
  
  base <- .build_tree_base(hc, title = NULL)
  p_tree <- base$plot
  tips   <- base$tips
  leaf_ids <- hc$labels
  
  # --- inner ring: limma beta
  beta_df <- limma_st %>%
    filter(protein %in% leaf_ids) %>%
    transmute(label = protein, beta = as.numeric(beta))
  
  bmax <- max(abs(beta_df$beta), na.rm = TRUE)
  if (!is.finite(bmax) || bmax == 0) bmax <- 0.1
  
  ring0 <- max(tips$x, na.rm = TRUE)
  thick <- 0.08
  gap   <- 0.02
  r1    <- ring0 + 0.05
  r2    <- r1 + thick + gap
  r_lab <- r2 + thick + 0.10
  
  p_tree <- .add_ring(
    p_tree, tips, beta_df,
    value_col  = "beta",
    ring_inner = r1, ring_thick = thick,
    fill_scale = scale_fill_gradient2(
      low = "#2166AC", mid = "grey97", high = "#B2182B",
      midpoint = 0,
      limits = c(-bmax, bmax),
      oob = squish,
      name = "limma beta (per 10y)"
    ),
    legend_title = "limma beta (per 10y)",
    show_legend = FALSE
  )
  
  # --- outer ring: cluster assignment (your hard-coded sets)
  prot_short <- .strip_human(prot)
  prot_up    <- toupper(prot_short)
  
  if (sex_label == "male") {
    cluster1_set <- toupper(c("cd5l","lv469","ibp2","pafa","pigr","saa4","lg3bp","ldha","vasp"))
    cluster <- ifelse(prot_up %in% cluster1_set, "1", "2")
  } else {
    cluster1_set <- toupper(c("ibp2","inhbc","lysc","igj","prop"))
    cluster <- ifelse(prot_up %in% cluster1_set, "1", "2")
  }
  
  cluster_df <- tibble(
    label   = prot,
    cluster = factor(cluster, levels = c("1","2"))
  )
  
  cluster_cols <- get_cluster_cols(sex_label)
  
  p_tree <- .add_ring(
    p_tree, tips, cluster_df,
    value_col  = "cluster",
    ring_inner = r2, ring_thick = thick,
    fill_scale = scale_fill_manual(
      values = cluster_cols,
      drop   = FALSE,
      name   = "Cluster"
    ),
    legend_title = "Cluster",
    show_legend = FALSE
  )
  
  # tip labels
  p_tree <- .add_tip_labels(p_tree, tips, r_label = r_lab, label_size = label_size)
  
  # widen canvas
  p_tree + xlim(0, r_lab + 0.15)
}

legend_cluster_shades <- function(cluster_cols,
                                  title = "Cluster",
                                  key_size = 0.18,
                                  title_size = 8,
                                  text_size  = 7) {
  
  k <- length(cluster_cols)
  xs <- 1 + (seq_len(k) - 1) * 2.1
  
  df <- data.frame(
    key  = factor(names(cluster_cols), levels = names(cluster_cols)),
    x    = xs,
    y    = 1.0,
    fill = unname(cluster_cols)
  )
  
  ggplot() +
    geom_tile(
      data = df,
      aes(x = x, y = y, fill = fill),
      width  = key_size, height = key_size,
      color  = "grey30", linewidth = 0.35,
      show.legend = FALSE
    ) +
    geom_text(
      data = df,
      aes(x = x, y = y, label = as.character(key)),
      hjust = 0,
      nudge_x = key_size/2 + 0.18,
      size = text_size / ggplot2::.pt
    ) +
    scale_fill_identity() +
    coord_fixed(
      ratio = 1,
      xlim  = c(0.55, max(xs) + 1.5),
      ylim  = c(0.65, 1.35),
      clip  = "off"
    ) +
    theme_void(base_size = 8, base_family = PLOT_FONT) +
    theme(
      plot.title  = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(title = title)
}

legend_limma_beta <- function(title = "limma beta (per 10 years)",
                              x_min = -0.9, x_max = 0.9,
                              bar_height = 0.045,
                              n = 360,
                              title_size = 8,
                              axis_text  = 7) {
  
  df <- data.frame(
    x = seq(x_min, x_max, length.out = n),
    y = 1,
    z = seq(x_min, x_max, length.out = n)
  )
  
  ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile(height = bar_height) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey97", high = "#B2182B",
      midpoint = 0,
      limits = c(x_min, x_max),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = c(x_min, 0, x_max),
      labels = c(sprintf("%.1f", x_min), "0", sprintf("%.1f", x_max)),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(0.92, 1.08), clip = "off") +
    theme_classic(base_size = 8, base_family = PLOT_FONT) +
    theme(
      plot.title   = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      axis.title   = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.line.x  = element_blank(),
      axis.text.x  = element_text(family = PLOT_FONT, size = axis_text),
      axis.ticks.x = element_line(color = "grey20"),
      panel.grid   = element_blank(),
      plot.margin  = margin(0, 0, 0, 0)
    ) +
    labs(title = title)
}

force_closed_circle_by_y_limits <- function(p) {
  gb <- ggplot_build(p)
  
  collect_y <- function(d) {
    out <- numeric(0)
    if (!is.null(d) && nrow(d) > 0) {
      if ("y"    %in% names(d)) out <- c(out, d$y)
      if ("ymin" %in% names(d)) out <- c(out, d$ymin)
      if ("ymax" %in% names(d)) out <- c(out, d$ymax)
    }
    out
  }
  
  ys <- unlist(lapply(gb$data, collect_y))
  ys <- ys[is.finite(ys)]
  if (length(ys) < 2) return(p)
  
  yr <- range(ys, na.rm = TRUE)
  p + scale_y_continuous(limits = yr, expand = expansion(mult = 0))
}

# ---- build circular panels (NO tags inside) ----
p_circ_m <- build_circular_panel_2rings("male",   hc_method = "average", label_size = 2.5) +
  labs(title = "Male") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = 10,
                              margin = margin(t = 3, r = 0, b = 3, l = 5)),
    plot.title.position = "plot"
  )

p_circ_f <- build_circular_panel_2rings("female", hc_method = "average", label_size = 2.5) +
  labs(title = "Female") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = 10,
                              margin = margin(t = 3, r = 0, b = 3, l = 5)),
    plot.title.position = "plot"
  )

p_circ_m <- force_closed_circle_by_y_limits(p_circ_m)
p_circ_f <- force_closed_circle_by_y_limits(p_circ_f)

# legends
leg_clust_m <- legend_cluster_shades(get_cluster_cols("male"),   title = "Cluster")
leg_clust_f <- legend_cluster_shades(get_cluster_cols("female"), title = "Cluster")
leg_beta    <- legend_limma_beta(title = "limma beta (per 10 years)")

row_top   <- (p_circ_m + p_circ_f) + plot_layout(ncol = 2, widths = c(1, 1)) &
  theme(plot.margin = margin(3, 5, 2, 5))

row_clust <- (leg_clust_m + leg_clust_f) + plot_layout(ncol = 2, widths = c(1, 1)) &
  theme(plot.margin = margin(0, 5, 0, 5))

row_beta  <- (plot_spacer() + leg_beta + plot_spacer()) +
  plot_layout(ncol = 3, widths = c(1.25, 1.00, 1.25)) &
  theme(plot.margin = margin(0, 5, 3, 5))

p_circ_block <- row_top / row_clust / row_beta +
  plot_layout(heights = c(1.00, 0.12, 0.07))

## ============================================================
## 3) FINAL COMBINATION (A/B at BLOCK LEVEL ONLY)
## ============================================================
p_combined_final <-
  (wrap_elements(full = p_conf_block) /
     wrap_elements(full = p_circ_block)) +
  plot_layout(heights = c(1.00, 1.15)) +
  plot_annotation(tag_levels = "A") &   # -> A for top block, B for bottom block
  theme(
    plot.tag = element_text(family = PLOT_FONT, size = 12, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(3, 5, 3, 5)
  )

p_combined_final

save_main_gg(
  p_combined_final,
  file.path(fig_main_dir, "confusion_circle"),
  height_mm = 190
)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 0) Stable proteins (>= 2/3)
# ------------------------------------------------------------
stopifnot(exists("res_male_nested"), exists("res_female_nested"))

stable_m <- res_male_nested$stability %>%
  filter(selection_freq >= 2/3) %>%
  pull(protein) %>% unique() %>%
  sub("_HUMAN$", "", .) %>%
  toupper() %>% sort()

stable_f <- res_female_nested$stability %>%
  filter(selection_freq >= 2/3) %>%
  pull(protein) %>% unique() %>%
  sub("_HUMAN$", "", .) %>%
  toupper() %>% sort()

force_common <- c("LDHA", "IBP2")
common_prots <- intersect(intersect(stable_m, stable_f), force_common)

male_only    <- setdiff(stable_m, common_prots)
female_only  <- setdiff(stable_f, common_prots)

# ------------------------------------------------------------
# 1) Geometry helpers
# ------------------------------------------------------------
circle_df <- function(cx, cy, r, n = 720) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
}
dist2 <- function(x, y, cx, cy) (x - cx)^2 + (y - cy)^2

# ------------------------------------------------------------
# 2) TEXT-SAFE margins
# ------------------------------------------------------------
text_margin <- function(label,
                        base = 0.22,      # baseline inward margin
                        per_char = 0.010, # extra per character
                        cap = 0.34) {     # cap per non esagerare
  m <- base + per_char * nchar(label)
  pmin(m, cap)
}

inside_with_label_margin <- function(x, y, cx, cy, r, label) {
  m <- text_margin(label)
  dist2(x, y, cx, cy) <= (r - m)^2
}

outside_with_label_margin <- function(x, y, cx, cy, r, label) {
  m <- text_margin(label)
  dist2(x, y, cx, cy) >= (r + m)^2
}

# ------------------------------------------------------------
# 3) 
# ------------------------------------------------------------
place_labels_textsafe <- function(labels, region,
                                  cx_m, cx_f, cy, r,
                                  xlim, ylim,
                                  grid_nx = 90, grid_ny = 64,
                                  seed = 1) {
  
  labels <- as.character(labels)
  n <- length(labels)
  if (n == 0) {
    return(tibble(label = character(0), x = numeric(0), y = numeric(0)))
  }
  
  gx <- seq(xlim[1], xlim[2], length.out = grid_nx)
  gy <- seq(ylim[1], ylim[2], length.out = grid_ny)
  G  <- expand.grid(x = gx, y = gy)
  
  set.seed(seed)
  used <- data.frame(x = numeric(0), y = numeric(0))
  out  <- vector("list", n)
  
  labels_ord <- labels
  
  for (i in seq_len(n)) {
    lab <- labels_ord[i]
    
    in_m  <- inside_with_label_margin(G$x, G$y, cx_m, cy, r, lab)
    in_f  <- inside_with_label_margin(G$x, G$y, cx_f, cy, r, lab)
    out_m <- outside_with_label_margin(G$x, G$y, cx_m, cy, r, lab)
    out_f <- outside_with_label_margin(G$x, G$y, cx_f, cy, r, lab)
    
    ok <- switch(
      region,
      male_only   = in_m & out_f,  
      female_only = in_f & out_m,  
      common      = in_m & in_f    
    )
    
    cand <- G[ok, , drop = FALSE]
    
    if (nrow(cand) < 5) {
      m2 <- function(label) {
        base <- 0.20; per_char <- 0.008; cap <- 0.30
        pmin(base + per_char * nchar(label), cap)
      }
      mm <- m2(lab)
      
      in_m2  <- dist2(G$x, G$y, cx_m, cy) <= (r - mm)^2
      in_f2  <- dist2(G$x, G$y, cx_f, cy) <= (r - mm)^2
      out_m2 <- dist2(G$x, G$y, cx_m, cy) >= (r + mm)^2
      out_f2 <- dist2(G$x, G$y, cx_f, cy) >= (r + mm)^2
      
      ok <- switch(
        region,
        male_only   = in_m2 & out_f2,
        female_only = in_f2 & out_m2,
        common      = in_m2 & in_f2
      )
      cand <- G[ok, , drop = FALSE]
    }
    
    if (nrow(cand) == 0) {
      xy <- switch(
        region,
        male_only   = c(cx_m - 0.55, cy),
        female_only = c(cx_f + 0.55, cy),
        common      = c(0,           cy)
      )
      out[[i]] <- data.frame(label = lab, x = xy[1], y = xy[2])
      next
    }
    
    if (nrow(used) > 0) {
      score <- apply(cand, 1, function(row) {
        dx <- row["x"] - used$x
        dy <- row["y"] - used$y
        min(dx * dx + dy * dy)
      })
      cand <- cand[order(score, decreasing = TRUE), , drop = FALSE]
    } else {
      cx0 <- mean(cand$x); cy0 <- mean(cand$y)
      d0  <- (cand$x - cx0)^2 + (cand$y - cy0)^2
      cand <- cand[order(d0), , drop = FALSE]
    }
    
    pick <- cand[1, , drop = FALSE]
    used <- rbind(used, pick)
    out[[i]] <- data.frame(label = lab, x = pick$x, y = pick$y)
  }
  
  bind_rows(out)
}

# ------------------------------------------------------------
# 4) Layout + circles
# ------------------------------------------------------------
male_col   <- "#009E73"
female_col <- "#E69F00"

cx_m <- -0.85
cx_f <-  0.85
cy   <-  0.12
r    <-  1.25

circ_m <- circle_df(cx_m, cy, r)
circ_f <- circle_df(cx_f, cy, r)

bounds <- list(
  male_only   = list(xlim = c(-2.35,  0.05), ylim = c(-1.10, 1.65)),
  female_only = list(xlim = c(-0.05,  2.35), ylim = c(-1.10, 1.65))
)

male_only_ord   <- sort(pub_protein(male_only))
female_only_ord <- sort(pub_protein(female_only))
common_ord      <- sort(pub_protein(common_prots))  # only LDHA/IBP2

df_m <- place_labels_textsafe(
  male_only_ord, "male_only",
  cx_m, cx_f, cy, r,
  bounds$male_only$xlim, bounds$male_only$ylim,
  seed = 101
)

df_f <- place_labels_textsafe(
  female_only_ord, "female_only",
  cx_m, cx_f, cy, r,
  bounds$female_only$xlim, bounds$female_only$ylim,
  seed = 303
)

df_common <- tibble(
  region = "Shared",
  label  = common_ord,
  x      = 0,
  y      = (cy + c(0.20, -0.20))[seq_along(common_ord)]
)

df_lab <- bind_rows(
  df_m %>% mutate(region = "Male only"),
  df_common,
  df_f %>% mutate(region = "Female only")
)

nice_family <- PLOT_FONT
n_total     <- nrow(df_lab)
txt_size    <- if (n_total <= 22) 5.6 else if (n_total <= 32) 5.0 else 4.4

# ------------------------------------------------------------
# 5) Plot
# ------------------------------------------------------------
p_venn <- ggplot() +
  geom_polygon(
    data = circ_m, aes(x = x, y = y),
    fill = male_col, alpha = 0.22,
    color = male_col, linewidth = 1.2
  ) +
  geom_polygon(
    data = circ_f, aes(x = x, y = y),
    fill = female_col, alpha = 0.22,
    color = female_col, linewidth = 1.2
  ) +
  geom_text(
    data = df_lab,
    aes(x = x, y = y, label = label),
    family = nice_family,
    size   = txt_size
  ) +
  annotate(
    "text", x = cx_m, y = -1.35, label = "Male",
    fontface = "bold", size = 5.4, family = nice_family
  ) +
  annotate(
    "text", x = cx_f, y = -1.35, label = "Female",
    fontface = "bold", size = 5.4, family = nice_family
  ) +
  coord_equal(xlim = c(-2.45, 2.45), ylim = c(-1.60, 1.85), clip = "off") +
  theme_void(base_size = 12, base_family = PLOT_FONT) +
  theme(plot.margin = margin(10, 18, 10, 18))

save_gg(
  p_venn,
  file.path(fig_supp_dir, "venn_diagram"),
  width_px = 1374, height_px = 862
)

# Figure S4: full protein-protein correlation heatmaps with stable-protein annotations.
save_heatmap_plot <- function(draw_fun, stem,
                              width_px, height_px, dpi = 96) {
  
  png_file <- paste0(stem, ".png")
  pdf_file <- paste0(stem, ".pdf")
  
  # Cartella temporanea LOCALE, non sul drive G:
  tmp_dir <- tempdir()
  
  png_tmp <- tempfile(pattern = "heatmap_", tmpdir = tmp_dir, fileext = ".png")
  
  message("PNG temporaneo: ", png_tmp)
  
  # Same ragg PNG device and local temporary-file sequence as the working version.
  ragg::agg_png(
    filename = png_tmp,
    width = width_px,
    height = height_px,
    units = "px",
    res = dpi,
    background = "white"
  )
  png_device_id <- grDevices::dev.cur()
  
  tryCatch(
    {
      draw_fun()
    },
    error = function(e) {
      message("Errore nel disegno PNG: ", conditionMessage(e))
      stop(e)
    },
    finally = {
      if (png_device_id %in% grDevices::dev.list()) {
        grDevices::dev.off(which = png_device_id)
      }
    }
  )
  
  # Controllo dopo la chiusura del device
  if (!file.exists(png_tmp)) {
    stop(
      "Il device PNG non ha creato il file temporaneo. Percorso atteso: ",
      png_tmp
    )
  }
  
  png_size <- file.info(png_tmp)$size
  
  if (is.na(png_size) || png_size <= 0) {
    stop("Il file PNG temporaneo è vuoto: ", png_tmp)
  }
  
  message("PNG temporaneo creato: ", png_size, " bytes")
  
  # Copia dal disco locale al drive G:
  copied <- file.copy(
    from = png_tmp,
    to = png_file,
    overwrite = TRUE
  )
  
  if (!copied || !file.exists(png_file)) {
    stop(
      "Il PNG è stato creato localmente, ma non può essere copiato nella destinazione: ",
      png_file
    )
  }
  
  # PDF: Cairo resolves Arial; draw/close sequence as before.
  grDevices::cairo_pdf(
    filename = pdf_file,
    width = width_px / dpi,
    height = height_px / dpi,
    family = PLOT_FONT,
    bg = "white"
  )
  pdf_device_id <- grDevices::dev.cur()
  
  tryCatch(
    {
      draw_fun()
    },
    finally = {
      if (pdf_device_id %in% grDevices::dev.list()) {
        grDevices::dev.off(which = pdf_device_id)
      }
    }
  )
  
  if (!file.exists(pdf_file) || is.na(file.info(pdf_file)$size) ||
      file.info(pdf_file)$size <= 0) {
    stop("Il file PDF non è stato completato: ", pdf_file)
  }
  
  message("PNG finale: ", normalizePath(png_file, winslash = "/"))
  message("Dimensione PNG finale: ", file.info(png_file)$size, " bytes")
  message("PDF finale: ", normalizePath(pdf_file, winslash = "/"))
  
  invisible(c(
    png = png_file,
    pdf = pdf_file
  ))
}

# Apply Arial only AFTER pheatmap has completed its internal
# PDF(NULL) layout pass. Passing the system font directly to pheatmap can fail
# in that temporary standard-PDF device before anything is drawn to the output.
.heatmap_times_new_roman <- function(g) {
  if (inherits(g, "text")) {
    if (is.null(g$gp)) g$gp <- grid::gpar()
    g$gp$fontfamily <- PLOT_FONT
  }
  if (length(g$grobs)) {
    g$grobs <- lapply(g$grobs, .heatmap_times_new_roman)
  }
  if (length(g$children)) {
    g$children <- do.call(grid::gList,
                          lapply(g$children, .heatmap_times_new_roman))
  }
  g
}

make_cor_heatmap <- function(sex_label,
                             stable_tbl,
                             stable_col,
                             selected_only = FALSE) {
  samples <- pheno_min$sample_id[pheno_min$sex == sex_label]
  samples <- intersect(samples, colnames(expr_clean))
  
  if (selected_only) {
    protein_ids <- stable_tbl$protein[stable_tbl$selection_freq >= 2/3]
    protein_ids <- intersect(protein_ids, rownames(expr_clean))
    stopifnot(length(protein_ids) >= 2)
    X <- expr_clean[protein_ids, samples, drop = FALSE]
  } else {
    X <- expr_clean[, samples, drop = FALSE]
  }
  
  C <- suppressWarnings(cor(t(X), use = "pairwise.complete.obs", method = "pearson"))
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  
  # Publication labels are applied directly to the selected correlation matrix.
  # This avoids a separate labels_row/labels_col rendering path.
  if (selected_only) {
    lab <- pub_protein(rownames(C))
    rownames(C) <- lab
    colnames(C) <- lab
  }
  
  if (!selected_only) {
    stable_ids <- stable_tbl$protein[stable_tbl$selection_freq >= 2/3]
    ann <- data.frame(
      Stable = factor(ifelse(rownames(C) %in% stable_ids, "Stable", "Other"),
                      levels = c("Other", "Stable")),
      row.names = rownames(C)
    )
    ann_cols <- list(Stable = c("Other" = "#F2F2F2", "Stable" = stable_col))
  } else {
    ann <- NA
    ann_cols <- NA
  }
  
  heatmap_grob <- pheatmap::pheatmap(
    C,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = selected_only,
    show_colnames = selected_only,
    annotation_row = ann,
    annotation_col = ann,
    annotation_colors = ann_cols,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    breaks = seq(-1, 1, length.out = 102),
    border_color = NA,
    main = ifelse(sex_label == "male", "Male", "Female"),
    fontsize = 10,
    angle_col = 45,
    silent = TRUE
  )$gtable
  
  .heatmap_times_new_roman(heatmap_grob)
}

# Full-proteome heatmap: same output as before.
draw_heatmap_pair <- function() {
  output_device <- grDevices::dev.cur()
  gm <- make_cor_heatmap(
    "male", stable_prots_m, "#009E73",
    selected_only = FALSE
  )
  
  gf <- make_cor_heatmap(
    "female", stable_prots_f, "#E69F00",
    selected_only = FALSE
  )
  
  # pheatmap's internal PDF(NULL) device can change the active device.
  # Explicitly return to the PNG/PDF opened by save_heatmap_plot().
  if (!(output_device %in% grDevices::dev.list())) {
    stop("The heatmap output device was closed unexpectedly.")
  }
  grDevices::dev.set(which = output_device)
  
  grid::grid.newpage()
  
  lay <- grid::grid.layout(1, 2)
  vp <- grid::viewport(layout = lay)
  
  grid::pushViewport(vp)
  
  grid::grid.draw(gridExtra::arrangeGrob(
    gm, gf,
    ncol = 2,
    top = NULL
  ))
  
  grid::grid.text(
    "A",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    gp = grid::gpar(fontfamily = PLOT_FONT, fontsize = 22, fontface = "bold")
  )
  
  grid::grid.text(
    "B",
    x = grid::unit(0.515, "npc"),
    y = grid::unit(0.965, "npc"),
    gp = grid::gpar(fontfamily = PLOT_FONT, fontsize = 22, fontface = "bold")
  )
  
  grid::popViewport()
}

save_heatmap_plot(
  draw_heatmap_pair,
  file.path(fig_supp_dir, "heatmap"),
  width_px = 1374, height_px = 862
)

# Additional supplementary figure restricted to the proteins selected by
# elastic net (>= 2/3 outer folds). No Stable/Other annotation bars.
draw_selected_heatmap_pair <- function() {
  output_device <- grDevices::dev.cur()
  gm <- make_cor_heatmap(
    "male", stable_prots_m, "#009E73",
    selected_only = TRUE
  )
  
  gf <- make_cor_heatmap(
    "female", stable_prots_f, "#E69F00",
    selected_only = TRUE
  )
  
  # pheatmap's internal PDF(NULL) device can change the active device.
  # Explicitly return to the PNG/PDF opened by save_heatmap_plot().
  if (!(output_device %in% grDevices::dev.list())) {
    stop("The heatmap output device was closed unexpectedly.")
  }
  grDevices::dev.set(which = output_device)
  
  grid::grid.newpage()
  
  lay <- grid::grid.layout(1, 2)
  vp <- grid::viewport(layout = lay)
  
  grid::pushViewport(vp)
  
  grid::grid.draw(gridExtra::arrangeGrob(
    gm, gf,
    ncol = 2,
    top = NULL
  ))
  
  grid::grid.text(
    "A",
    x = grid::unit(0.025, "npc"),
    y = grid::unit(0.965, "npc"),
    gp = grid::gpar(fontfamily = PLOT_FONT, fontsize = 22, fontface = "bold")
  )
  
  grid::grid.text(
    "B",
    x = grid::unit(0.515, "npc"),
    y = grid::unit(0.965, "npc"),
    gp = grid::gpar(fontfamily = PLOT_FONT, fontsize = 22, fontface = "bold")
  )
  
  grid::popViewport()
}

selected_heatmap_stem <- file.path(fig_supp_dir, "heatmap_selected_proteins")

save_heatmap_plot(
  draw_selected_heatmap_pair,
  selected_heatmap_stem,
  width_px = 1374, height_px = 862
)

# Fail explicitly if the PNG was not completed successfully.
selected_png <- paste0(selected_heatmap_stem, ".png")
if (!file.exists(selected_png) || is.na(file.info(selected_png)$size) ||
    file.info(selected_png)$size <= 0) {
  stop("heatmap_selected_proteins.png was not written correctly. Check the R console for the graphics-device error immediately above this message.")
}

## ============================================================
## Line plots (STANDARDIZED): mean across age classes (Male vs Female)
##   - Panel A: IBP2
##   - Panel B: LDHA
## Standardization: per-protein z-score across ALL samples (male + female)
## Variability: mean ± SE ribbons within each sex x ageclass
##
## Requires: expr_clean (proteins x samples), pheno_min (aligned)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---- 0) Defensive checks ----
stopifnot(is.matrix(expr_clean))
stopifnot(all(c("sample_id","sex","ageclass") %in% names(pheno_min)))
stopifnot(all(pheno_min$sample_id %in% colnames(expr_clean)))

# ---- 1) Robustly find protein IDs in expr_clean ----
find_prot_id <- function(mat_rownames, prot) {
  hit <- grep(paste0("^", prot, "(_HUMAN)?$"), mat_rownames, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) stop("Protein not found in expr_clean rownames: ", prot)
  hit[1]
}

ibp2_id <- find_prot_id(rownames(expr_clean), "IBP2")
ldha_id <- find_prot_id(rownames(expr_clean), "LDHA")

# ---- 2) Extract named vectors + standardize (z-score) ----
extract_named_row <- function(mat, row_id) {
  v <- mat[row_id, , drop = TRUE]          # numeric vector
  v <- as.vector(v)                        # keep as plain vector
  names(v) <- colnames(mat)                # force names = sample IDs (always)
  v
}

standardize_z <- function(v) {
  # v must be numeric; names are preserved
  nm <- names(v)
  mu <- mean(v, na.rm = TRUE)
  sdv <- stats::sd(v, na.rm = TRUE)
  
  if (!is.finite(sdv) || sdv <= 0) {
    out <- rep(0, length(v))
    names(out) <- nm
    return(out)
  }
  
  out <- (v - mu) / sdv
  names(out) <- nm
  out
}

x_ibp2 <- extract_named_row(expr_clean, ibp2_id)
x_ldha <- extract_named_row(expr_clean, ldha_id)

Z_ibp2 <- standardize_z(x_ibp2)
Z_ldha <- standardize_z(x_ldha)

# ---- 3) Summaries: mean ± SE by sex x ageclass ----
make_mean_se_df_from_vec <- function(z_vec, prot_label) {
  # Align to pheno_min order
  z_aligned <- z_vec[pheno_min$sample_id]
  
  df <- tibble(
    sample_id = pheno_min$sample_id,
    sex       = pheno_min$sex,
    ageclass  = pheno_min$ageclass,
    z         = as.numeric(z_aligned)
  ) %>%
    filter(is.finite(z))
  
  if (nrow(df) == 0) stop("No finite values after alignment for protein: ", prot_label)
  
  df %>%
    group_by(sex, ageclass) %>%
    summarise(
      mean_z = mean(z),
      sd_z   = sd(z),
      n      = n(),
      se_z   = sd_z / sqrt(n),
      ymin   = mean_z - se_z,
      ymax   = mean_z + se_z,
      .groups = "drop"
    ) %>%
    mutate(protein = prot_label)
}

df_ibp2 <- make_mean_se_df_from_vec(Z_ibp2, "IBP2")
df_ldha <- make_mean_se_df_from_vec(Z_ldha, "LDHA")

# ---- 4) Plot helper ----
SEX_COL <- c(male = "#009E73", female = "#E69F00")

plot_mean_lines_with_var_std <- function(df, title, show_legend = TRUE,
                                         base_size = 12,
                                         axis_text = 12.5,
                                         axis_title = 13.5,
                                         title_size = 18) {
  
  ggplot(df, aes(x = ageclass, y = mean_z, color = sex, fill = sex, group = sex)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.8) +
    scale_x_discrete(
      labels = c(
        "20-29" = "20-29",
        "30-39" = "30-39",
        "40-49" = "40-49",
        "50-59" = "50-59",
        "60-70" = "60-69"
      )
    ) +
    scale_color_manual(values = SEX_COL, name = NULL) +
    scale_fill_manual(values = SEX_COL, name = NULL) +
    labs(
      title = title,
      x = "Age class",
      y = "Mean intensity (cleaned proteins, standardized)"
    ) +
    theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    theme(
      plot.title   = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      axis.text.x  = element_text(family = PLOT_FONT, size = axis_text),
      axis.text.y  = element_text(family = PLOT_FONT, size = axis_text),
      axis.title.x = element_text(family = PLOT_FONT, size = axis_title),
      axis.title.y = element_text(family = PLOT_FONT, size = axis_title),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text  = element_text(family = PLOT_FONT, size = axis_text),
      plot.margin  = margin(10, 25, 10, 25)
    )
}

p_ibp2 <- plot_mean_lines_with_var_std(df_ibp2, pub_protein("IBP2"), show_legend = TRUE) + labs(tag = "A")
p_ldha <- plot_mean_lines_with_var_std(df_ldha, "LDHA", show_legend = TRUE) + labs(tag = "B")

# ---- 5) Combine into A/B panels with shared legend ----
p_two_prots_std <-
  (p_ibp2 + p_ldha) +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(
    plot.tag = element_text(family = PLOT_FONT, size = 22, face = "bold"),
    legend.position = "bottom"
  )

p_two_prots_std

save_gg(
  p_two_prots_std,
  file.path(fig_supp_dir, "ibp2_ldha_plot"),
  width_px = 1060, height_px = 862
)

# 7. AgeProt, Figure S6 and descriptive tables
# --- paths

expr_path <- here("analysis","data","proteomic_cleaned_expr.rds")
meta_path <- here("analysis","data","proteomic_nako_metadata_aligned.rds")
nako_path <- here("analysis","data","nako_base.rds")

stable_m_path <- here("analysis","results","protein_sets","stable_proteins_male_ids.rds")
stable_f_path <- here("analysis","results","protein_sets","stable_proteins_female_ids.rds")

# --- load

expr_clean <- as.matrix(readRDS(expr_path)) # proteins x samples
meta_raw <- readRDS(meta_path) # aligned metadata
nako_base <- readRDS(nako_path) # contains Sample_ID + a_rr_hypertens

stable_ids_m <- readRDS(stable_m_path) # character vector
stable_ids_f <- readRDS(stable_f_path)

# --- helper: clean numeric variables with typical NAKO missing codes

clean_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  bad_codes <- c(-99, -98, -97, -88, -9, 7775, 7777, 9999)
  x[x %in% bad_codes] <- NA_real_
  x
}

# --- pick age-in-years column from meta_raw (robust)

age_col <- intersect(c("age_years","age","Age","AGE"), names(meta_raw))[1]
if (is.na(age_col)) {
  stop("Nessuna colonna di età in anni trovata in meta_raw. Attese: age_years / age / Age / AGE.")
}

# 1) Minimal pheno from aligned metadata (include age in years)

pheno <- meta_raw %>%
  transmute(
    sample_id    = as.character(Sample_ID),
    sex          = factor(sex, levels = c(1,2), labels = c("male","female")),
    age_years    = clean_num(.data[[age_col]]),
    age10        = age_years / 10,  # age in decades (effect per 10 years)
    ageclass     = droplevels(as.factor(ageclass)),  # kept for reference/plots if needed
    bmi          = as.numeric(bmi),
    smok_status  = fct_explicit_na(as.factor(smok_status),  na_level = "Missing"),
    edu          = fct_explicit_na(as.factor(edu),          na_level = "Missing"),
    alcohol_cons = fct_explicit_na(as.factor(alcohol_cons), na_level = "Missing")
  ) %>%
  mutate(
    bmi = ifelse(is.na(bmi), median(bmi, na.rm = TRUE), bmi),
    # reference levels (modifica se vuoi)
    smok_status  = fct_relevel(smok_status, "non-smoker"),
    edu          = fct_relevel(edu, "medium"),
    alcohol_cons = fct_relevel(alcohol_cons, "none")
  )

# 2) Add selected NAKO covariates from nako_base (via Sample_ID)

nako_cov <- as.data.frame(nako_base) %>%
  transmute(
    sample_id        = as.character(Sample_ID),
    a_rr_hypertens   = clean_num(a_rr_hypertens),
    a_rr_sys         = clean_num(a_rr_sys),
    a_rr_dia         = clean_num(a_rr_dia),
    a_gpaq_ptotalday = clean_num(a_gpaq_ptotalday),
    weist            = clean_num(anthro_taillenumfang)#,
    #diabetes         = as.factor(d_an_met_1)
  )

pheno <- pheno %>%
  left_join(nako_cov, by = "sample_id") %>%
  mutate(
    hypertens = case_when(
      a_rr_hypertens == 0 ~ "No",
      a_rr_hypertens == 1 ~ "Yes",
      TRUE                ~ "Missing"
    ),
    hypertens = factor(hypertens, levels = c("No","Yes","Missing"))
  )

# 3) Align pheno to expression columns

pheno <- pheno %>%
  filter(sample_id %in% colnames(expr_clean)) %>%
  arrange(match(sample_id, colnames(expr_clean)))

stopifnot(identical(pheno$sample_id, colnames(expr_clean)))

# ============================================================



# ============================================================

covars <- c(
  "age10"#,                 # age in decades (chronological age)
  # "bmi_z",
  # "smok_status",
  # 'a_rr_sys',
  # 'a_rr_dia',
  # 'a_gpaq_ptotalday',
  # 'weist'#,
  #'diabetes'
)

compute_ageprot_score <- function(expr_mat, pheno_df, stable_ids, sex_label) {
  
  ## samples for this sex (aligned to expr columns)
  sel_samples <- pheno_df$sample_id[pheno_df$sex == sex_label]
  sel_samples <- intersect(sel_samples, colnames(expr_mat))
  
  stable_prots <- intersect(stable_ids, rownames(expr_mat))
  
  stopifnot(length(sel_samples) >= 3)
  stopifnot(length(stable_prots) >= 2)
  
  ## expression: proteins x samples
  X <- expr_mat[stable_prots, sel_samples, drop = FALSE]
  
  ## pheno aligned to X columns
  ph <- pheno_df %>%
    dplyr::filter(sex == sex_label, sample_id %in% sel_samples) %>%
    dplyr::arrange(match(sample_id, sel_samples))
  
  stopifnot(identical(ph$sample_id, sel_samples))
  
  design <- model.matrix(~ age10, data = ph)
  rownames(design) <- ph$sample_id
  stopifnot(identical(colnames(X), rownames(design)))
  
  ## 1) standardize each protein across subjects (z-score per protein)
  X_z <- t(scale(t(X), center = TRUE, scale = TRUE))
  
  ## 2) estimate coefficients on the same (z) scale
  fit <- limma::eBayes(limma::lmFit(X_z, design))
  beta_age <- fit$coefficients[, "age10"]
  
  ## guard against NA (e.g., constant proteins)
  beta_age[is.na(beta_age)] <- 0
  
  ## 3) compute score on the same standardized scale
  AgeProt_score <- as.numeric(t(X_z) %*% beta_age)
  names(AgeProt_score) <- colnames(X_z)
  
  ## return score per sample (aligned to sel_samples)
  tibble::tibble(sample_id = sel_samples, AgeProt_score = AgeProt_score)
}

run_ageprot_lm <- function(sex_label, stable_ids, covars, expr_mat, pheno_df) {
  
  scores <- compute_ageprot_score(expr_mat, pheno_df, stable_ids, sex_label)
  
  df <- pheno_df %>%
    dplyr::filter(sex == sex_label) %>%
    dplyr::inner_join(scores, by = "sample_id") %>%
    dplyr::mutate(
      ## within-sex standardizations (continuous)
      bmi_z              = as.numeric(scale(bmi)),
      a_rr_sys_z         = as.numeric(scale(a_rr_sys)),
      a_rr_dia_z         = as.numeric(scale(a_rr_dia)),
      a_gpaq_ptotalday_z = as.numeric(scale(a_gpaq_ptotalday)),
      ## drop unused levels (categorical)
      smok_status  = droplevels(smok_status),
      edu          = droplevels(edu),
      alcohol_cons = droplevels(alcohol_cons),
      hypertens    = droplevels(hypertens)
    )
  
  ## check covariates exist
  missing_vars <- setdiff(covars, names(df))
  if (length(missing_vars) > 0) {
    stop(
      paste0(
        "Covariate mancanti in df: ", paste(missing_vars, collapse = ", "),
        "\nDisponibili: ", paste(names(df), collapse = ", ")
      )
    )
  }
  
  ## keep complete cases on outcome + covariates
  df <- df %>%
    dplyr::filter(complete.cases(dplyr::select(., AgeProt_score, dplyr::all_of(covars))))
  
  fit <- stats::lm(stats::reformulate(covars, response = "AgeProt_score"), data = df)
  
  list(data = df, fit = fit)
}

## --- male
res_male <- run_ageprot_lm(
  sex_label = "male",
  stable_ids = stable_ids_m,
  covars = covars,
  expr_mat = expr_clean,
  pheno_df = pheno
)

cat("\n========================\nMALE: model summary\n========================\n")
print(summary(res_male$fit))

## --- female
res_female <- run_ageprot_lm(
  sex_label = "female",
  stable_ids = stable_ids_f,
  covars = covars,
  expr_mat = expr_clean,
  pheno_df = pheno
)

cat("\n========================\nFEMALE: model summary\n========================\n")
print(summary(res_female$fit))

## ============================================================
## Scatter: AgeProt_score ~ age10 (separate panels by sex)
## - outcome: AgeProt_score
## - predictor: age10 (chronological age, decades)
## - descriptive / internal validation
## ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

COL_M <- "#009E73"
COL_F <- "#E69F00"

## lighter hues for points + ribbons
PT_M  <- "#66C2A5"
PT_F  <- "#F6B26B"
RIB_M <- "#BFE7D9"
RIB_F <- "#FBE0B8"

plot_ageprot_scatter_panel <- function(df, sex_label,
                                       line_col, point_col, ribbon_col,
                                       axis_text  = 12.5,
                                       axis_title = 13.5,
                                       title_size = 18,
                                       base_size  = 12) {
  
  d <- df %>%
    filter(sex == sex_label) %>%
    filter(is.finite(AgeProt_score), is.finite(age10))
  
  stopifnot(nrow(d) >= 3)
  
  ## AgeProt as OUTCOME
  fit <- lm(AgeProt_score ~ age10, data = d)
  r2  <- summary(fit)$r.squared
  r2_txt <- sprintf("R² = %.2f", r2)
  
  ## annotation position (top-left)
  x_rng <- range(d$age10, na.rm = TRUE)
  y_rng <- range(d$AgeProt_score, na.rm = TRUE)
  x_pos <- x_rng[1] + 0.05 * diff(x_rng)
  y_pos <- y_rng[2] - 0.07 * diff(y_rng)
  
  ggplot(d, aes(x = age10, y = AgeProt_score)) +
    geom_point(size = 2.6, alpha = 0.45, color = point_col) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 1.05,
      color = line_col,
      fill  = ribbon_col
    ) +
    annotate(
      "text",
      x = x_pos, y = y_pos,
      label = r2_txt,
      hjust = 0, vjust = 1,
      fontface = "bold",
      size = 4.6
    ) +
    labs(
      title = ifelse(sex_label == "male", "Male", "Female"),
      x = "Age (decades)",
      y = "AgeProt score"
    ) +
    theme_classic(base_size = base_size, base_family = PLOT_FONT) +
    theme(
      plot.title   = element_text(family = PLOT_FONT, face = "bold", hjust = 0, size = title_size),
      axis.text.x  = element_text(family = PLOT_FONT, size = axis_text),
      axis.text.y  = element_text(family = PLOT_FONT, size = axis_text),
      axis.title.x = element_text(family = PLOT_FONT, size = axis_title),
      axis.title.y = element_text(family = PLOT_FONT, size = axis_title),
      plot.margin  = margin(10, 25, 10, 25)
    )
}

## combined data
df_plot <- bind_rows(res_male$data, res_female$data)

## panels
p_m <- plot_ageprot_scatter_panel(
  df_plot, "male",
  line_col = COL_M, point_col = PT_M, ribbon_col = RIB_M
) + labs(tag = "A")

p_f <- plot_ageprot_scatter_panel(
  df_plot, "female",
  line_col = COL_F, point_col = PT_F, ribbon_col = RIB_F
) + labs(tag = "B")

## final layout
p_ageprot_scatter <-
  (p_m + p_f) +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.margin   = margin(10, 25, 10, 25),
    panel.spacing = unit(2.0, "cm"),
    plot.tag      = element_text(family = PLOT_FONT, size = 22, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

p_ageprot_scatter

save_gg(
  p_ageprot_scatter,
  file.path(fig_supp_dir, "ageprotscore_scatter"),
  width_px = 992, height_px = 758
)

library(dplyr)
library(forcats)
library(tidyr)
library(gt)

# ------------------------------------------------------------
# Goal: avoid creating a tiny 70–79 class due to a single 70.22.
# Strategy: cap ages at 69.99 before building decade labels.
# This forces 70.22 -> 69.99 -> "60-69".
# ------------------------------------------------------------

ph_tab <- pheno %>%
  mutate(
    age_years = suppressWarnings(as.numeric(age_years)),
    age_years_cap = dplyr::if_else(!is.na(age_years) & age_years >= 70, 69.99, age_years),
    
    age_dec = floor(age_years_cap / 10) * 10,
    ageclass_lbl = ifelse(is.na(age_dec), NA_character_, paste0(age_dec, "-", age_dec + 9)),
    
    bmi = suppressWarnings(as.numeric(bmi)),
    bmi_cat = case_when(
      is.na(bmi) ~ "Missing",
      bmi < 18.5 ~ "Underweight",
      bmi < 25   ~ "Normal",
      bmi < 30   ~ "Overweight",
      TRUE       ~ "Obese"
    ) %>% factor(levels = c("Underweight", "Normal", "Overweight", "Obese", "Missing")),
    
    sex          = droplevels(sex),
    smok_status  = fct_explicit_na(as.factor(smok_status),  na_level = "Missing") %>% droplevels(),
    alcohol_cons = fct_explicit_na(as.factor(alcohol_cons), na_level = "Missing") %>% droplevels(),
    edu          = fct_explicit_na(as.factor(edu),          na_level = "Missing") %>% droplevels()
  )

# Define *fixed* age levels you want to show (no 70-79)
age_levels <- c("20-29","30-39","40-49","50-59","60-69")
ph_tab <- ph_tab %>%
  mutate(ageclass_lbl = factor(ageclass_lbl, levels = age_levels))

safe_colname <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
age_levels_safe <- safe_colname(age_levels)

# ----------------------------
# Helpers: n (%) by ageclass within sex + Overall (no p-values)
# ----------------------------

fmt_np <- function(n, p) sprintf("%d (%.1f%%)", n, p)

summ_cat_by <- function(df, var, by) {
  tab <- table(df[[var]], df[[by]], useNA = "ifany")
  pct <- prop.table(tab, 2) * 100
  
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  out$pct  <- as.data.frame(pct, stringsAsFactors = FALSE)$Freq
  out$cell <- mapply(fmt_np, out$Freq, out$pct)
  
  names(out)[1:2] <- c("level", "group")
  out[, c("level", "group", "cell")]
}

summ_cat_overall <- function(df, var) {
  tab <- table(df[[var]], useNA = "ifany")
  pct <- prop.table(tab) * 100
  tibble(
    level = names(tab),
    overall = mapply(fmt_np, as.integer(tab), as.numeric(pct))
  )
}

build_sex_block <- function(df_sex, sex_prefix, vars, by = "ageclass_lbl") {
  df_sex <- df_sex %>%
    mutate(ageclass_lbl = factor(.data[[by]], levels = age_levels)) %>%
    droplevels()
  
  pieces <- lapply(vars, function(v) {
    long <- summ_cat_by(df_sex, v, by) %>%
      pivot_wider(names_from = group, values_from = cell)
    
    overall <- summ_cat_overall(df_sex, v)
    
    out <- long %>%
      left_join(overall, by = "level") %>%
      mutate(
        Variable = v,
        Level = as.character(level)
      )
    
    names(out) <- safe_colname(names(out))
    
    for (nm in age_levels_safe) {
      if (!nm %in% names(out)) out[[nm]] <- NA_character_
    }
    
    out %>%
      select(Variable, Level, all_of(age_levels_safe), overall) %>%
      rename_with(~ paste0(sex_prefix, "_", .x), -c(Variable, Level))
  })
  
  bind_rows(pieces)
}

# ----------------------------
# Build one single GT table with column spanners (Male / Female)
# ----------------------------

vars <- c("bmi_cat", "smok_status", "alcohol_cons", "edu")

df_m <- ph_tab %>% filter(sex == "male") %>% droplevels()
df_f <- ph_tab %>% filter(sex == "female") %>% droplevels()

male_df   <- build_sex_block(df_m, "Male", vars, by = "ageclass_lbl")
female_df <- build_sex_block(df_f, "Female", vars, by = "ageclass_lbl")

tab_df <- male_df %>%
  full_join(female_df, by = c("Variable", "Level")) %>%
  mutate(
    Variable = recode(
      Variable,
      bmi_cat       = "BMI category",
      smok_status   = "Smoking status",
      alcohol_cons  = "Alcohol consumption",
      edu           = "Education"
    )
  )

male_cols   <- c(paste0("Male_", age_levels_safe), "Male_overall")
female_cols <- c(paste0("Female_", age_levels_safe), "Female_overall")

male_labels <- setNames(c(age_levels, "Overall"), male_cols)
female_labels <- setNames(c(age_levels, "Overall"), female_cols)

gt_table <- tab_df %>%
  gt(groupname_col = "Variable") %>%
  tab_header(title = "Table 1. Participant characteristics by age class, stratified by sex") %>%
  cols_label(Level = "", !!!male_labels, !!!female_labels) %>%
  tab_spanner(label = "Male",   columns = all_of(male_cols)) %>%
  tab_spanner(label = "Female", columns = all_of(female_cols)) %>%
  sub_missing(everything(), missing_text = "") %>%
  cols_align(align = "left",   columns = c(Level)) %>%
  cols_align(align = "center", columns = c(all_of(male_cols), all_of(female_cols))) %>%
  tab_options(
    table.width = pct(100),
    table.font.size = px(10),
    row_group.font.weight = "bold",
    data_row.padding = px(2)
  ) %>%
  tab_footnote(
    footnote = "Values are n (%).",
    locations = cells_title(groups = "title")
  )

gt_table

# Optional save:
# out_dir <- here::here("analysis", "results", "tables")
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# gt::gtsave(gt_table, filename = file.path(out_dir, "Table1_by_ageclass_stratified_sex.html"))
# gt::gtsave(gt_table, filename = file.path(out_dir, "Table1_by_ageclass_stratified_sex.png"), expand = 0)
# try(gt::gtsave(gt_table, filename = file.path(out_dir, "Table1_by_ageclass_stratified_sex.rtf")))

save_gt(
  gt_table,
  file.path(tables_dir, "Table1_by_ageclass_stratified_sex"),
  vwidth = 1952, vheight = 960
)

# Main Table 1 (descriptive table used in the current manuscript code).
# Use observed metadata for the same analytical participants; do not change pheno,
# model covariates, imputations, or any preceding/supplementary output.

t1_required <- c("Sample_ID", "bmi", "edu", "smok_status",
                 "anthro_taillenumfang", "a_rr_sys", "a_rr_dia",
                 "a_rr_hypertens", "a_gpaq_ptotalday",
                 paste0("d_an_cv_", 1:6), "d_an_met_1", "d_an_ca")
if (!all(t1_required %in% names(nako_base))) {
  stop("Table 1: missing metadata columns: ",
       paste(setdiff(t1_required, names(nako_base)), collapse = ", "))
}
stopifnot(!anyNA(pheno$sample_id), !anyDuplicated(pheno$sample_id),
          !anyNA(nako_base$Sample_ID), !anyDuplicated(nako_base$Sample_ID))
t1_idx <- match(pheno$sample_id, as.character(nako_base$Sample_ID))
stopifnot(!anyNA(t1_idx))
t1_raw <- as.data.frame(nako_base)[t1_idx, , drop = FALSE]
stopifnot(identical(as.character(t1_raw$Sample_ID), pheno$sample_id))

t1_num <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[!is.finite(x) | x %in% c(-99, -98, -97, -88, -9, -5,
                             7775, 7777, 8888, 9999,
                             777777, 888888, 999999)] <- NA_real_
  x
}
t1_category <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) %in% c("", "Missing", "NA")] <- NA_character_
  x
}
t1_yes_no <- function(x, yes = 1, no = 2) {
  # Anamnesis: 1 = yes, 2 = no; unknown/not specified stays missing.
  x <- t1_num(x)
  ifelse(x %in% c(yes, no), x == yes, NA)
}

t1_cvd <- vapply(t1_raw[paste0("d_an_cv_", 1:5)],
                 t1_yes_no, logical(nrow(t1_raw)))
# Any positive component establishes CVD; absence requires all five known negatives.
t1_cvd_any <- ifelse(rowSums(t1_cvd, na.rm = TRUE) > 0, TRUE,
                     ifelse(rowSums(!is.na(t1_cvd)) == ncol(t1_cvd), FALSE, NA))

ph_main <- tibble(
  sample_id = pheno$sample_id,
  age_years = pheno$age_years,
  # Same group membership as before; >=60 also correctly labels ages >=70.
  age_group = cut(pheno$age_years, breaks = c(20, 30, 40, 50, 60, Inf),
                  right = FALSE, labels = c("20-29", "30-39", "40-49", "50-59", "\u226560")),
  female = as.character(pheno$sex) == "female",
  high_edu = t1_category(t1_raw$edu) == "high",
  non_smoker = t1_category(t1_raw$smok_status) == "non-smoker",
  former_smoker = t1_category(t1_raw$smok_status) == "former smoker",
  smoker = t1_category(t1_raw$smok_status) == "smoker",
  bmi = t1_num(t1_raw$bmi),
  bmi_gt30 = bmi > 30,
  waist = t1_num(t1_raw$anthro_taillenumfang),
  systolic_bp = t1_num(t1_raw$a_rr_sys),
  diastolic_bp = t1_num(t1_raw$a_rr_dia),
  measured_hypertension = t1_yes_no(t1_raw$a_rr_hypertens, yes = 1, no = 0),
  activity = t1_num(t1_raw$a_gpaq_ptotalday),
  reported_hypertension = t1_yes_no(t1_raw$d_an_cv_6),
  cardiovascular_disease = t1_cvd_any,
  diabetes = t1_yes_no(t1_raw$d_an_met_1),
  cancer = t1_yes_no(t1_raw$d_an_ca)
)
stopifnot(nrow(ph_main) == nrow(pheno), !anyNA(ph_main$age_group))

groups <- c("Total", levels(ph_main$age_group))
t1_spec <- tibble::tribble(
  ~variable, ~Characteristic, ~summary_type,
  "age_years", "Age, years, mean (SD)", "mean_sd",
  "female", "Female sex, n (%)", "n_pct",
  "high_edu", "High educational attainment, n (%)", "n_pct",
  "non_smoker", "Non-smoker, n (%)", "n_pct",
  "former_smoker", "Former smoker, n (%)", "n_pct",
  "smoker", "Smoker, n (%)", "n_pct",
  "bmi", "BMI (kg/m²), mean (SD)", "mean_sd",
  "bmi_gt30", "BMI >30 kg/m², n (%)", "n_pct",
  "waist", "Waist circumference, cm, mean (SD)", "mean_sd",
  "systolic_bp", "Systolic blood pressure, mmHg, mean (SD)", "mean_sd",
  "diastolic_bp", "Diastolic blood pressure, mmHg, mean (SD)", "mean_sd",
  "measured_hypertension", "Elevated measured blood pressure, n (%)", "n_pct",
  "activity", "Physical activity (GPAQ), min/day, median [Q1, Q3]", "median_iqr",
  "reported_hypertension", "History of diagnosed hypertension, n (%)", "n_pct",
  "cardiovascular_disease", "History of cardiovascular disease, n (%)", "n_pct",
  "diabetes", "History of diabetes mellitus, n (%)", "n_pct",
  "cancer", "History of cancer, n (%)", "n_pct"
)
t1_format <- function(x, summary_type) {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  switch(summary_type,
         mean_sd = if (length(x) > 1L) sprintf("%.1f (%.1f)", mean(x), sd(x))
         else sprintf("%.1f (NA)", mean(x)),
         median_iqr = sprintf("%.1f [%.1f, %.1f]", median(x),
                              quantile(x, 0.25, names = FALSE),
                              quantile(x, 0.75, names = FALSE)),
         n_pct = sprintf("%d (%.1f%%)", sum(x), 100 * mean(x))
  )
}
table1_df <- t1_spec %>% select(Characteristic)
table1_df$Available <- vapply(t1_spec$variable,
                              function(v) sum(!is.na(ph_main[[v]])), integer(1))
for (g in groups) {
  t1_subset <- if (g == "Total") ph_main else filter(ph_main, age_group == g)
  table1_df[[g]] <- vapply(seq_len(nrow(t1_spec)), function(i) {
    t1_format(t1_subset[[t1_spec$variable[i]]], t1_spec$summary_type[i])
  }, character(1))
}

n_by_group <- c(Total = nrow(ph_main), table(ph_main$age_group))
label_map <- setNames(
  paste0(groups, "\n(N=", as.integer(n_by_group[groups]), ")"), groups
)
gt_table1_main <- table1_df %>%
  gt() %>%
  tab_header(title = "Table 1. Characteristics of the study population by age group") %>%
  cols_label(Characteristic = "Characteristics", Available = "Available, n", !!!label_map) %>%
  cols_align("left", Characteristic) %>%
  cols_align("center", all_of(c("Available", groups))) %>%
  tab_options(table.width = pct(100), table.font.size = px(10), data_row.padding = px(2)) %>%
  tab_source_note(source_note = paste(
    "Values are mean (SD), median [Q1, Q3], or n (%).",
    "Available, n denotes non-missing observations in the total sample;",
    "percentages use the non-missing denominator within each column.",
    "Unknown/not specified responses and NAKO missing codes are excluded; no imputation is used.",
    "High educational attainment corresponds to the 'high' metadata category.",
    "BMI, body mass index; SD, standard deviation; Q1/Q3, 25th/75th percentiles."
  )) %>%
  tab_source_note(source_note = paste(
    "Elevated measured blood pressure: systolic >=140 or diastolic >=90 mmHg",
    "(NAKO a_rr_hypertens), without considering antihypertensive medication.",
    "Disease histories are self-reported physician diagnoses; cardiovascular disease includes",
    "myocardial infarction, coronary artery disease/angina, heart failure, arrhythmia,",
    "or peripheral arterial disease, with hypertension reported separately."
  )) %>%
  tab_source_note(source_note = paste(
    "GPAQ, Global Physical Activity Questionnaire; total physical activity across work,",
    "transport and leisure, in minutes/day. The oldest age group includes all participants aged >=60 years."
  ))

save_gt(
  gt_table1_main,
  file.path(tables_dir, "Table1_main"),
  vwidth = 1952, vheight = 1100
)

## ============================================================
## Helpers
## ============================================================

# 8. Figure 4: partial-correlation networks
prepare_sex_data_stable <- function(sex_label,
                                    stable_prots,
                                    expr_mat,
                                    pheno_df,
                                    min_samples  = 20,
                                    min_proteins = 5,
                                    scale_z      = TRUE) {
  stopifnot(sex_label %in% levels(pheno_df$sex))
  expr_mat <- as.matrix(expr_mat)
  
  sel_samples <- pheno_df$sample_id[pheno_df$sex == sex_label]
  sel_samples <- intersect(sel_samples, colnames(expr_mat))
  
  if (length(sel_samples) < min_samples) {
    stop("Too few samples for sex='", sex_label, "': n=", length(sel_samples))
  }
  
  stable_in_expr <- intersect(unique(stable_prots), rownames(expr_mat))
  if (length(stable_in_expr) < min_proteins) {
    stop("Too few stable proteins found in expr for sex='", sex_label,
         "': p=", length(stable_in_expr))
  }
  
  ## bootnet wants rows=samples, cols=variables
  X <- t(expr_mat[stable_in_expr, sel_samples, drop = FALSE])
  
  ## drop zero-variance proteins
  v <- apply(X, 2, var, na.rm = TRUE)
  keep <- is.finite(v) & v > 0
  X <- X[, keep, drop = FALSE]
  
  if (ncol(X) < min_proteins) {
    stop("After variance filter, too few proteins remain for sex='", sex_label,
         "': p=", ncol(X))
  }
  
  if (scale_z) X <- scale(X)
  
  X <- as.data.frame(X)
  cat("\n[", sex_label, "] data ready: n=", nrow(X), " p=", ncol(X), "\n", sep = "")
  return(X)
}

is_empty_graph <- function(net_obj) {
  W <- net_obj$graph
  if (is.null(W)) return(TRUE)
  diag(W) <- 0
  all(W == 0)
}

## ============================================================
## Pipeline: pcor + nonparam bootstrap + bootThreshold
##  - default="pcor" evita gamma/lambda tuning (niente EBICglasso)
##  - alpha=1 e threshold=0: NON filtro a priori => rete non vuota
##  - bootThreshold(boots, alpha=0.05): rendo sparsa via CI bootstrap
## ============================================================
run_pcor_bootThreshold <- function(dat,
                                   corMethod = "cor",    # o "cor_auto"/"npn"
                                   nCores = max(1, parallel::detectCores() - 1),
                                   nBoots_nonparam = 500,
                                   thr_alpha = 0.05,     # livello per bootThreshold
                                   computeCentrality = FALSE) {
  
  net <- estimateNetwork(
    dat,
    default   = "pcor",
    corMethod = corMethod,
    alpha     = 1,
    threshold = 0
  )
  
  ##  nonparam bootstrap for CI edge
  set.seed(12345)
  boots <- bootnet(
    net,
    nBoots  = nBoots_nonparam,
    type    = "nonparametric",
    nCores  = 1,
    statistics = c("edge"),
    computeCentrality = computeCentrality
  )
  
  net_thr <- bootThreshold(boots, alpha = thr_alpha)
  
  list(
    data = dat,
    net = net,
    boot_nonparam = boots,
    net_thresholded = net_thr
  )
}

edge_table_from_net <- function(net_obj) {
  W <- net_obj$graph
  diag(W) <- 0
  idx <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return(tibble(from = character(), to = character(), weight = numeric()))
  }
  tibble(
    from = rownames(W)[idx[, 1]],
    to   = colnames(W)[idx[, 2]],
    weight = W[idx]
  ) %>% arrange(desc(abs(weight)))
}

## ============================================================
## Run male / female (aging-driven)
## ============================================================

CORMETHOD <- "cor"
NCORES <- max(1, parallel::detectCores() - 1)

NBOOT_NONPARAM <- 500
THR_ALPHA <- 0.1

## --- male ---
dat_male <- prepare_sex_data_stable(
  sex_label    = "male",
  stable_prots = stable_ids_m,
  expr_mat     = expr_clean,
  pheno_df     = pheno_min,
  min_samples  = 30,
  scale_z      = TRUE
)

net_res_male <- run_pcor_bootThreshold(
  dat_male,
  corMethod = CORMETHOD,
  nCores = NCORES,
  nBoots_nonparam = NBOOT_NONPARAM,
  thr_alpha = THR_ALPHA,
  computeCentrality = FALSE
)

saveRDS(net_res_male, file = file.path(network_dir, "aging_driven_pcor_male.rds"))

## --- female ---
dat_female <- prepare_sex_data_stable(
  sex_label    = "female",
  stable_prots = stable_ids_f,
  expr_mat     = expr_clean,
  pheno_df     = pheno_min,
  min_samples  = 30,
  scale_z      = TRUE
)

net_res_female <- run_pcor_bootThreshold(
  dat_female,
  corMethod = CORMETHOD,
  nCores = NCORES,
  nBoots_nonparam = NBOOT_NONPARAM,
  thr_alpha = THR_ALPHA,
  computeCentrality = FALSE
)

saveRDS(net_res_female, file = file.path(network_dir, "aging_driven_pcor_female.rds"))

cat("\nSaved network results to:\n", network_dir, "\n", sep = "")

library(igraph)

SEX_COL <- c(male = "#009E73", female = "#E69F00")

strip_human_suffix <- function(x) {
  x <- as.character(x)
  x <- gsub("(_HUMAN|\\|HUMAN|\\sHUMAN)$", "", x, ignore.case = TRUE, perl = TRUE)
  x <- gsub("[_\\|\\s]+$", "", x, perl = TRUE)
  trimws(x)
}

get_weight_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x) && !is.null(x$graph)) {
    g <- x$graph
    if (is.matrix(g)) return(g)
    if (is.data.frame(g)) return(as.matrix(g))
  }
  if (is.list(x) && !is.null(x$results) && is.matrix(x$results)) return(x$results)
  if (is.list(x) && !is.null(x$network)  && is.matrix(x$network))  return(x$network)
  stop("Cannot extract a weight matrix from object of class: ", paste(class(x), collapse = ", "))
}

sanitize_W <- function(W) {
  W <- as.matrix(W)
  suppressWarnings(storage.mode(W) <- "numeric")
  if (nrow(W) != ncol(W)) stop("Weight matrix is not square: ", nrow(W), "x", ncol(W))
  
  if (is.null(colnames(W)) && !is.null(rownames(W))) colnames(W) <- rownames(W)
  if (is.null(rownames(W)) && !is.null(colnames(W))) rownames(W) <- colnames(W)
  if (is.null(rownames(W)) && is.null(colnames(W))) {
    nm <- paste0("V", seq_len(ncol(W)))
    rownames(W) <- nm
    colnames(W) <- nm
  }
  
  diag(W) <- 0
  W[!is.finite(W)] <- 0
  W
}

extract_W_and_labels <- function(net_obj) {
  W <- sanitize_W(get_weight_matrix(net_obj))
  
  labs <- NULL
  if (inherits(net_obj, "qgraph") && !is.null(net_obj$labels)) labs <- net_obj$labels
  if (is.null(labs) && is.list(net_obj) && !is.null(net_obj$labels)) labs <- net_obj$labels
  
  if (is.null(labs)) labs <- colnames(W)
  
  if (length(labs) != ncol(W)) labs <- colnames(W)
  
  list(W = W, labels = labs)
}

# Common edge-strength scale across male and female panels.
# Edge SIGN is encoded only by colour (red = positive, blue = negative).
# Edge MAGNITUDE is encoded only by line width (thicker = larger |partial correlation|).
get_nonzero_abs_weights <- function(net_obj) {
  W <- extract_W_and_labels(net_obj)$W
  v <- abs(W[upper.tri(W) & W != 0])
  v[is.finite(v)]
}

all_edge_abs <- c(
  get_nonzero_abs_weights(net_res_male$net_thresholded),
  get_nonzero_abs_weights(net_res_female$net_thresholded)
)

if (length(all_edge_abs) == 0) {
  NETWORK_EDGE_LIMITS <- c(0, 1)
} else {
  NETWORK_EDGE_LIMITS <- range(all_edge_abs, na.rm = TRUE)
  if (!all(is.finite(NETWORK_EDGE_LIMITS))) NETWORK_EDGE_LIMITS <- c(0, 1)
  if (diff(NETWORK_EDGE_LIMITS) < 1e-12) {
    NETWORK_EDGE_LIMITS <- c(0, NETWORK_EDGE_LIMITS[2])
  }
}

scale_edge_strength <- function(absw, limits = NETWORK_EDGE_LIMITS) {
  if (length(absw) == 0) return(numeric(0))
  if (!all(is.finite(limits)) || diff(limits) <= 0) return(rep(1, length(absw)))
  pmin(1, pmax(0, (absw - limits[1]) / diff(limits)))
}

edge_width_by_strength <- function(weights,
                                   limits = NETWORK_EDGE_LIMITS,
                                   min_edge_width = 1.0,
                                   max_edge_width = 7.0) {
  s <- scale_edge_strength(abs(weights), limits)
  min_edge_width + s * (max_edge_width - min_edge_width)
}

plot_thresholded_network_simple <- function(net_obj,
                                            sex_label = c("male", "female"),
                                            panel_letter = "A",
                                            title_text = NULL,
                                            vertex_size = 30,
                                            label_cex = 0.85,
                                            frame_width = 1.2,
                                            min_edge_width = 1.0,
                                            max_edge_width = 7.0,
                                            edge_abs_limits = NETWORK_EDGE_LIMITS,
                                            layout_seed = 123,
                                            margin = 0.25) {
  
  sex_label <- match.arg(sex_label)
  col_frame <- SEX_COL[sex_label]
  
  obj <- extract_W_and_labels(net_obj)
  W <- obj$W
  labs <- pub_protein(strip_human_suffix(obj$labels))
  
  if (all(W == 0)) {
    plot.new()
    mtext(panel_letter, side = 3, line = 1.2, adj = 0.06, font = 2, cex = 1.5)
    if (is.null(title_text)) title_text <- tools::toTitleCase(sex_label)
    mtext(title_text, side = 3, line = 0.2, adj = 0.35, font = 2, cex = 1.25)
    text(0.5, 0.5, "Empty after bootThreshold", cex = 1)
    return(invisible(NULL))
  }
  
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  idx <- match(igraph::V(g)$name, colnames(W))
  igraph::V(g)$label <- labs[idx]
  node_label_cex <- rep(label_cex, igraph::vcount(g))
  if (sex_label == "female") {
    node_label_cex[igraph::V(g)$label == "IGFBP2"] <- 0.65#label_cex * (0.75 / 0.85)
  }
  
  w <- igraph::E(g)$weight
  absw <- abs(w)
  
  centrality <- data.frame(
    protein = igraph::V(g)$label,
    degree = igraph::degree(g),
    strength = igraph::strength(g, weights = absw),
    betweenness = igraph::betweenness(g, directed = FALSE, weights = 1 / absw)
  )
  
  print(centrality[order(-centrality$strength), ])
  
  # Sign and magnitude have deliberately separate visual encodings:
  # red/blue = sign; line width = absolute partial correlation.
  igraph::E(g)$color <- ifelse(w < 0, "#2166AC", "#B2182B")
  igraph::E(g)$width <- edge_width_by_strength(
    w,
    limits = edge_abs_limits,
    min_edge_width = min_edge_width,
    max_edge_width = max_edge_width
  )
  
  set.seed(layout_seed)
  lay <- igraph::layout_with_fr(g, weights = absw + 1e-6)
  
  if (sex_label == "female") {
    i <- match("IGFBP2", igraph::V(g)$label)
    
    if (!is.na(i)) {
      xr <- range(lay[, 1])
      yr <- range(lay[, 2])
      
      lay[i, 1] <- min(lay[i, 1] + 0.15 * diff(xr), xr[2])
      lay[i, 2] <- min(lay[i, 2] + 0.08 * diff(yr), yr[2])
    }
  }
  
  if (sex_label == "female") {
    i <- match("CP", igraph::V(g)$label)
    
    if (!is.na(i)) {
      xr <- range(lay[, 1])
      yr <- range(lay[, 2])
      
      lay[i, 1] <- min(lay[i, 1] + 0.04 * diff(xr), xr[2])
      lay[i, 2] <- max(lay[i, 2] - 0.04 * diff(yr), yr[1])
    }
  }
  
  if (sex_label == "female") {
    i <- match("GP9", igraph::V(g)$label)
    
    if (!is.na(i)) {
      xr <- range(lay[, 1])
      yr <- range(lay[, 2])
      
      lay[i, 1] <- min(lay[i, 1] + 0.04 * diff(xr), xr[2])
      lay[i, 2] <- max(lay[i, 2] - 0.04 * diff(yr), yr[1])
    }
  }
  
  op <- par(xpd = NA)
  on.exit(par(op), add = TRUE)
  
  plot(
    g,
    layout = lay,
    rescale = TRUE,
    margin = margin,
    vertex.color = "white",
    vertex.size = 0.93 * pmax(vertex_size, 5 * nchar(igraph::V(g)$label)),
    vertex.frame.color = col_frame,
    vertex.frame.width = frame_width,
    vertex.label = igraph::V(g)$label,
    vertex.label.cex = node_label_cex,
    vertex.label.family = PLOT_FONT,
    vertex.label.color = "black",
    vertex.label.dist = 0,
    edge.color = igraph::E(g)$color,
    edge.width = igraph::E(g)$width,
    edge.label = NA,
    edge.curved = 0.10,
    main = ""
  )
  
  mtext(panel_letter, side = 3, line = 1.2, adj = 0.06, font = 2, cex = 1.5)
  if (is.null(title_text)) title_text <- tools::toTitleCase(sex_label)
  mtext(title_text, side = 3, line = 0.2, adj = 0.35, font = 2, cex = 1.25)
  
  invisible(list(graph = g, W = W, layout = lay, labels_used = labs))
}

# Shared horizontal legend below the two network panels.
# Colour encodes the sign of the partial correlation; line width encodes its absolute magnitude.
draw_network_edge_legend <- function(edge_abs_limits = NETWORK_EDGE_LIMITS,
                                     min_edge_width = 1.0,
                                     max_edge_width = 7.0) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  lim <- edge_abs_limits
  if (!all(is.finite(lim)) || diff(lim) <= 0) lim <- c(0, 1)
  
  # One definition for the whole legend.
  
  # ----- sign key -----
  
  segments(0.03, 0.65, 0.09, 0.65, col = "#B2182B", lwd = 3.2, lend = "round")
  text(0.105, 0.65, "Positive", adj = c(0, 0.5), cex = 0.95)
  
  segments(0.23, 0.65, 0.29, 0.65, col = "#2166AC", lwd = 3.2, lend = "round")
  text(0.305, 0.65, "Negative", adj = c(0, 0.5), cex = 0.95)
  
  # ----- magnitude key: horizontal tapered bar -----
  text(0.52, 0.88, "Edge weight: |partial correlation|", adj = c(0, 0.5),
       font = 2, cex = 0.95)
  
  x0 <- 0.55
  x1 <- 0.97
  yc <- 0.55
  half_low  <- 0.010
  half_high <- 0.060
  
  polygon(
    x = c(x0, x1, x1, x0),
    y = c(yc - half_low, yc - half_high, yc + half_high, yc + half_low),
    col = "grey30", border = NA
  )
  
  tick_vals <- pretty(lim, n = 4)
  tick_vals <- tick_vals[tick_vals >= lim[1] - 1e-12 & tick_vals <= lim[2] + 1e-12]
  tick_vals <- sort(unique(c(lim[1], tick_vals, lim[2])))
  if (length(tick_vals) > 5) {
    tick_vals <- tick_vals[round(seq(1, length(tick_vals), length.out = 5))]
  }
  
  tick_x <- x0 + (tick_vals - lim[1]) / diff(lim) * (x1 - x0)
  for (i in seq_along(tick_vals)) {
    segments(tick_x[i], 0.43, tick_x[i], 0.47, col = "grey35", lwd = 0.8)
    text(tick_x[i], 0.34, sprintf("%.2f", tick_vals[i]),
         adj = c(0.5, 0.5), cex = 0.90)
  }
  
  invisible(NULL)
}

draw_network_figure <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  # Keep the two network panels at their original size and place the shared
  # legend in a separate horizontal strip below them.
  layout(
    matrix(c(1, 2,
             3, 3), nrow = 2, byrow = TRUE),
    widths = c(1, 1),
    heights = c(1, 0.17)
  )
  
  par(mar = c(1, 1, 3.2, 1), cex = 1)
  
  plot_thresholded_network_simple(
    net_res_male$net_thresholded,
    sex_label = "male",
    panel_letter = "A",
    title_text = "Male",
    vertex_size = 30,
    label_cex = 0.75,
    margin = 0.08,
    edge_abs_limits = NETWORK_EDGE_LIMITS
  )
  
  plot_thresholded_network_simple(
    net_res_female$net_thresholded,
    sex_label = "female",
    panel_letter = "B",
    title_text = "Female",
    vertex_size = 30,
    label_cex = 0.75,
    margin = 0.08,
    edge_abs_limits = NETWORK_EDGE_LIMITS
  )
  
  par(mar = c(0.3, 1, 0.2, 1), cex = 1)
  draw_network_edge_legend(NETWORK_EDGE_LIMITS)
}

save_main_base(
  draw_network_figure,
  file.path(fig_main_dir, "bootnet_results"),
  height_mm = 120
)

# 9. Figure 5: hallmark/pathway mapping
# ------------------------------------------------------------------
# Manual pathway assignment.
# IMPORTANT:
# - directions below are taken from the ordered barplot,
#   i.e. every protein has an explicit direction.
# - the figure is built as:
#     Official hallmark  ->  Protein marker
#   with full protein names written to the right of each protein symbol.
# ------------------------------------------------------------------

male_col   <- "#009E73"
female_col <- "#E69F00"
neutral_col <- "grey25"
PATHWAY_FONT <- PLOT_FONT

pathway_tbl <- tribble(
  ~sex,      ~protein, ~direction, ~pathway,                                       ~hallmark,
  "Male",   "IBP2",  "up",      "Insulin-like Growth Factor Signaling", "Deregulated Nutrient Sensing",
  "Male",   "LG3BP", "up",      "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Male",   "PIGR",  "up",      "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Male",   "VASP",  "up",      "Cytoskeletal Organization and Dynamics",     "Loss of Proteostasis",
  "Male",   "PAFA",  "up",      "Lipid Transport and Macrophage Signaling",   "Altered Intercellular Communication",
  "Male",   "LV469", "up",      "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Male",   "LDHA",  "up",      "Glycolysis and Glucose Metabolism",          "Deregulated Nutrient Sensing",
  "Male",   "SAA4",  "up",      "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Male",   "TENX",  "down",    "Extracellular Matrix Organization",    "Loss of Proteostasis",
  "Male",   "PI16",  "down",    "Extracellular Matrix Organization",    "Loss of Proteostasis",
  "Male",   "ATS13", "down",    "Extracellular Matrix Organization",    "Loss of Proteostasis",
  "Male",   "CD5L",  "down",    "Extracellular Matrix Organization",    "Loss of Proteostasis",
  "Male",   "IGA2",  "down",    "Humoral and Adaptive Immune Response",       "Chronic Inflammation",
  "Male",   "HV307", "down",    "Humoral and Adaptive Immune Response",       "Chronic Inflammation",
  "Male",   "CO6A1", "down",    "Extracellular Matrix Organization",    "Loss of Proteostasis",
  "Male",   "TRFL",  "down",    "Iron Homeostasis and Redox Regulation",      "Mitochondrial Dysfunction",
  "Male",   "LV211", "down",    "Humoral and Adaptive Immune Response",       "Chronic Inflammation",
  
  "Female", "IBP2",  "up",      "Insulin-like Growth Factor Signaling", "Deregulated Nutrient Sensing",
  "Female", "INHBC", "up",      "Insulin-like Growth Factor Signaling", "Deregulated Nutrient Sensing",
  "Female", "LYSC",  "up",      "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Female", "IGJ",   "up",      "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Female", "PROP",  "up",      "Insulin-like Growth Factor Signaling", "Deregulated Nutrient Sensing",
  "Female", "LDHA",  "down",    "Glycolysis and Glucose Metabolism",          "Deregulated Nutrient Sensing",
  "Female", "ALS",   "down",    "Insulin-like Growth Factor Signaling", "Deregulated Nutrient Sensing",
  "Female", "GPIX",  "down",    "Renin-Angiotensin and Vascular Regulation",  "Altered Intercellular Communication",
  "Female", "PRDX6", "down",    "Antioxidant Defense and Redox Homeostasis",  "Mitochondrial Dysfunction",
  "Female", "ITB3",  "down",    "Protease Inhibition and Chaperone Activity", "Loss of Proteostasis",
  "Female", "CERU",  "down",    "Antioxidant Defense and Redox Homeostasis",  "Mitochondrial Dysfunction",
  "Female", "ATPA",  "down",    "Antioxidant Defense and Redox Homeostasis",  "Mitochondrial Dysfunction",
  "Female", "HV315", "down",    "Acute Phase Response and Innate Immunity",   "Chronic Inflammation",
  "Female", "ANGT",  "down",    "Renin-Angiotensin and Vascular Regulation",  "Altered Intercellular Communication",
  "Female", "PZP",   "down",    "Protease Inhibition and Chaperone Activity", "Loss of Proteostasis"
)

hallmark_levels <- c(
  "Deregulated Nutrient Sensing",
  "Chronic Inflammation",
  "Loss of Proteostasis",
  "Altered Intercellular Communication",
  "Mitochondrial Dysfunction"
)


male_full_order <- c(
  "IBP2", "LG3BP", "PIGR", "VASP", "PAFA", "LV469", "LDHA", "SAA4",
  "TENX", "PI16", "ATS13", "CD5L", "IGA2", "HV307", "CO6A1", "TRFL", "LV211"
)

female_full_order <- c(
  "IBP2", "INHBC", "LYSC", "IGJ", "PROP", "LDHA", "ALS", "GPIX", "PRDX6",
  "ITB3", "CERU", "ATPA", "HV315", "ANGT", "PZP"
)

shared_proteins <- c("IBP2", "LDHA")

protein_levels <- c(
  setdiff(male_full_order, shared_proteins),
  shared_proteins,
  setdiff(female_full_order, shared_proteins)
)

arrow_from_direction <- function(x) {
  dplyr::case_when(
    x == "up"   ~ "↑",
    x == "down" ~ "↓",
    TRUE         ~ ""
  )
}

split_colored_protein <- function(protein, col_left, col_right) {
  split_idx <- floor(nchar(protein) / 2)
  left_part  <- substr(protein, 1, split_idx)
  right_part <- substr(protein, split_idx + 1, nchar(protein))
  
  paste0(
    "<span style='color:", col_left, "; font-weight:700;'>", left_part, "</span>",
    "<span style='color:", col_right, "; font-weight:700;'>", right_part, "</span>"
  )
}

plot_df <- pathway_tbl %>%
  mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    hallmark = factor(hallmark, levels = hallmark_levels),
    protein = factor(protein, levels = protein_levels)
  )

shared_dirs <- plot_df %>%
  filter(as.character(protein) %in% shared_proteins) %>%
  select(sex, protein, direction) %>%
  mutate(sex = as.character(sex)) %>%
  pivot_wider(names_from = sex, values_from = direction)

protein_meta <- plot_df %>%
  mutate(protein = as.character(protein)) %>%
  distinct(protein, pathway, hallmark) %>%
  left_join(shared_dirs, by = "protein") %>%
  mutate(
    group = case_when(
      protein %in% shared_proteins ~ "Shared",
      protein %in% setdiff(male_full_order, shared_proteins) ~ "Male",
      TRUE ~ "Female"
    ),
    pathway_wrap = stringr::str_wrap(pathway, width = 34),
    pathway_col = case_when(
      group == "Male"   ~ male_col,
      group == "Female" ~ female_col,
      TRUE               ~ neutral_col
    ),
    protein_label_html = case_when(
      group == "Male" ~ paste0(
        "<span style='color:", male_col, "; font-weight:700;'>", protein, "</span>",
        "&nbsp;&nbsp;",
        "<span style='color:", male_col, "; font-weight:700;'>", arrow_from_direction(Male), "</span>"
      ),
      group == "Female" ~ paste0(
        "<span style='color:", female_col, "; font-weight:700;'>", protein, "</span>",
        "&nbsp;&nbsp;",
        "<span style='color:", female_col, "; font-weight:700;'>", arrow_from_direction(Female), "</span>"
      ),
      protein == "IBP2" ~ paste0(
        split_colored_protein("IBP2", male_col, female_col),
        "&nbsp;&nbsp;",
        "<span style='color:", male_col, "; font-weight:700;'>", arrow_from_direction(Male), "</span>",
        "&nbsp;",
        "<span style='color:", female_col, "; font-weight:700;'>", arrow_from_direction(Female), "</span>"
      ),
      protein == "LDHA" ~ paste0(
        split_colored_protein("LDHA", male_col, female_col),
        "&nbsp;&nbsp;",
        "<span style='color:", male_col, "; font-weight:700;'>", arrow_from_direction(Male), "</span>",
        "&nbsp;",
        "<span style='color:", female_col, "; font-weight:700;'>", arrow_from_direction(Female), "</span>"
      ),
      TRUE ~ protein
    )
  )

# ============================================================
# TWO-PANEL SANKEY WITH CLASSIC RIBBONS
# - Male panel on top, Female panel below
# - panel titles in BLACK
# - protein names and pathway labels in BLACK
# - 5 bright hallmark colors (including green and orange)
# - hallmark labels equally spaced
# ============================================================

library(tidyverse)
library(patchwork)

# ---- 5 bright, distinct hallmark colors ----
hallmark_cols <- c(
  "Deregulated Nutrient Sensing" = "#2979FF",  # bright blue
  "Chronic Inflammation"         = "#FF1744",  # bright red
  "Loss of Proteostasis"         = "#D81BB2",  # bright violet
  "Altered Intercellular Communication"  = "#00C853",  # bright green
  "Mitochondrial Dysfunction"    = "#FF9100"   # bright orange
)

hallmark_levels <- c(
  "Deregulated Nutrient Sensing",
  "Chronic Inflammation",
  "Loss of Proteostasis",
  "Altered Intercellular Communication",
  "Mitochondrial Dysfunction"
)

# ---- order from the barplots ----
male_order_top <- c(
  "IBP2", "LG3BP", "PIGR", "VASP", "PAFA", "LV469", "LDHA", "SAA4",
  "TENX", "PI16", "ATS13", "CD5L", "IGA2", "HV307", "CO6A1", "TRFL", "LV211"
)

female_order_top <- c(
  "IBP2", "INHBC", "LYSC", "IGJ", "PROP", "LDHA", "ALS", "GPIX",
  "PRDX6", "ITB3", "CERU", "ATPA", "HV315", "ANGT", "PZP"
)

# Full protein names from the analysis protein annotation/manuscript.
# Keys retain the existing short UniProt IDs; displayed symbols, directions,
# hallmark assignments, row order and graphical/export settings are unchanged.
HALLMARK_PROTEIN_NAMES <- c(
  IBP2 = "Insulin-like growth factor-binding protein 2",
  LG3BP = "Galectin-3-binding protein",
  PIGR = "Polymeric immunoglobulin receptor",
  VASP = "Vasodilator-stimulated phosphoprotein",
  PAFA = "Platelet-activating factor acetylhydrolase",
  LV469 = "Immunoglobulin lambda variable 4-69",
  LDHA = "L-lactate dehydrogenase A chain",
  SAA4 = "Serum amyloid A-4 protein",
  TENX = "Tenascin-X",
  PI16 = "Peptidase inhibitor 16",
  ATS13 = "A disintegrin and metalloproteinase with thrombospondin motifs 13",
  CD5L = "CD5 antigen-like",
  IGA2 = "Immunoglobulin alpha-2 heavy chain",
  HV307 = "Immunoglobulin heavy variable 3-7",
  CO6A1 = "Collagen alpha-1(VI) chain",
  TRFL = "Lactotransferrin",
  LV211 = "Immunoglobulin lambda variable 2-11",
  INHBC = "Inhibin beta C chain",
  LYSC = "Lysozyme C",
  IGJ = "Immunoglobulin J chain",
  PROP = "Properdin",
  ALS = "Insulin-like growth factor-binding protein complex acid labile subunit",
  GPIX = "Platelet glycoprotein IX",
  PRDX6 = "Peroxiredoxin-6",
  ITB3 = "Integrin beta-3",
  CERU = "Ceruloplasmin",
  ATPA = "ATP synthase F(1) complex subunit alpha, mitochondrial",
  HV315 = "Immunoglobulin heavy variable 3-15",
  ANGT = "Angiotensinogen",
  PZP = "Pregnancy zone protein"
)
stopifnot(all(c(male_order_top, female_order_top) %in% names(HALLMARK_PROTEIN_NAMES)))

# ------------------------------------------------------------
# Ribbon helper
# ------------------------------------------------------------
bezier_curve <- function(p0, p1, p2, p3, t) {
  (1 - t)^3 * p0 +
    3 * (1 - t)^2 * t * p1 +
    3 * (1 - t) * t^2 * p2 +
    t^3 * p3
}

make_ribbon_polygon <- function(x0, x1, y0, y1, h = 0.24, n = 100) {
  t <- seq(0, 1, length.out = n)
  
  cx1 <- x0 + 0.33 * (x1 - x0)
  cx2 <- x0 + 0.67 * (x1 - x0)
  
  x_top <- bezier_curve(x0, cx1, cx2, x1, t)
  y_top <- bezier_curve(y0 + h/2, y0 + h/2, y1 + h/2, y1 + h/2, t)
  
  x_bot <- bezier_curve(x0, cx1, cx2, x1, t)
  y_bot <- bezier_curve(y0 - h/2, y0 - h/2, y1 - h/2, y1 - h/2, t)
  
  tibble(
    x = c(x_top, rev(x_bot)),
    y = c(y_top, rev(y_bot))
  )
}

# ------------------------------------------------------------
# Build one panel
# ------------------------------------------------------------
build_panel_objects <- function(sex_label, prot_order) {
  
  df <- pathway_tbl %>%
    filter(sex == sex_label, protein %in% prot_order) %>%
    mutate(
      protein = factor(protein, levels = prot_order),
      hallmark = factor(hallmark, levels = hallmark_levels),
      pathway_wrap = stringr::str_wrap(unname(HALLMARK_PROTEIN_NAMES[as.character(protein)]), width = 48),
      protein_label = paste0(pub_protein(protein), " ", arrow_from_direction(direction))
    ) %>%
    arrange(protein) %>%
    mutate(
      y_target = rev(seq_len(n()))
    )
  
  n_prot <- nrow(df)
  
  # equally spaced hallmark labels
  hallmark_df <- tibble(
    hallmark = factor(hallmark_levels, levels = hallmark_levels),
    y_label = seq(from = n_prot - 0.5, to = 1.5, length.out = length(hallmark_levels)),
    label = stringr::str_wrap(hallmark_levels, width = 25)
  )
  
  ribbon_h   <- 0.24
  ribbon_gap <- 0.045
  
  # source positions within each hallmark band
  df2 <- df %>%
    left_join(hallmark_df %>% select(hallmark, y_label), by = "hallmark") %>%
    group_by(hallmark) %>%
    arrange(desc(y_target), .by_group = TRUE) %>%
    mutate(
      n_h = n(),
      slot_index = row_number(),
      offset = ((n_h + 1) / 2 - slot_index) * (ribbon_h + ribbon_gap),
      y_source = y_label + offset
    ) %>%
    ungroup() %>%
    mutate(
      ribbon_h = ribbon_h
    )
  
  ribbon_df <- purrr::map2_dfr(
    seq_len(nrow(df2)),
    split(df2, seq_len(nrow(df2))),
    ~{
      row <- .y
      poly <- make_ribbon_polygon(
        x0 = 1.10, x1 = 1.86,
        y0 = row$y_source,
        y1 = row$y_target,
        h  = row$ribbon_h
      )
      poly %>%
        mutate(
          ribbon_id = .x,
          hallmark = row$hallmark
        )
    }
  )
  
  header_y <- max(df2$y_target) + 1.3
  
  list(
    row_df = df2,
    hallmark_df = hallmark_df,
    ribbon_df = ribbon_df,
    header_y = header_y
  )
}

male_obj   <- build_panel_objects("Male", male_order_top)
female_obj <- build_panel_objects("Female", female_order_top)

# ------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------
make_panel_plot <- function(obj, panel_title, show_legend = FALSE) {
  
  ggplot() +
    geom_polygon(
      data = obj$ribbon_df,
      aes(x = x, y = y, group = ribbon_id, fill = hallmark),
      color = NA,
      alpha = 0.88
    ) +
    geom_segment(
      data = obj$row_df,
      aes(x = 2.30, xend = 2.34, y = y_target, yend = y_target),
      linewidth = 0.18,
      color = "grey72"
    ) +
    geom_text(
      data = obj$hallmark_df,
      aes(x = 0.98, y = y_label, label = label, color = hallmark),
      hjust = 1,
      vjust = 0.5,
      fontface = "bold",
      family = PATHWAY_FONT,
      size = 7.5 / ggplot2::.pt,
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    geom_text(
      data = obj$row_df,
      aes(x = 2.06, y = y_target, label = protein_label),
      hjust = 0.5,
      vjust = 0.5,
      color = "black",
      family = PATHWAY_FONT,
      size = 7 / ggplot2::.pt
    ) +
    geom_text(
      data = obj$row_df,
      aes(x = 2.38, y = y_target, label = pathway_wrap),
      hjust = 0,
      vjust = 0.5,
      color = "black",
      family = PATHWAY_FONT,
      size = 7 / ggplot2::.pt,
      lineheight = 0.98
    ) +
    annotate(
      "text",
      x = 0.98, y = obj$header_y,
      label = "Hallmark of aging\n(López-Otín et al., 2023)",
      hjust = 1,
      fontface = "bold",
      family = PATHWAY_FONT,
      size = 7.5 / ggplot2::.pt
    ) +
    annotate(
      "text",
      x = 2.06, y = obj$header_y,
      label = "Protein\nmarker",
      fontface = "bold",
      family = PATHWAY_FONT,
      size = 7.5 / ggplot2::.pt
    ) +
    annotate(
      "text",
      x = 2.38, y = obj$header_y,
      label = "Protein name",
      hjust = 0,
      fontface = "bold",
      family = PATHWAY_FONT,
      size = 7.5 / ggplot2::.pt
    ) +
    scale_fill_manual(
      values = hallmark_cols,
      breaks = hallmark_levels,
      labels = hallmark_levels
    ) +
    scale_color_manual(
      values = hallmark_cols,
      breaks = hallmark_levels,
      labels = hallmark_levels
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, title = NULL,
                               keywidth = unit(3, "mm"),
                               keyheight = unit(3, "mm"))) +
    coord_cartesian(
      xlim = c(0.15, 3.60),
      ylim = c(0.5, obj$header_y + 1.0),
      clip = "off"
    ) +
    labs(
      title = panel_title,
      fill = NULL
    ) +
    theme_void(base_size = 8, base_family = PATHWAY_FONT) +
    theme(
      plot.title = element_text(
        face = "bold",
        family = PATHWAY_FONT,
        size = 10,
        color = "black",
        hjust = 0.12
      ),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 7, family = PATHWAY_FONT),
      legend.key.size = unit(3, "mm"),
      legend.spacing.x = unit(2, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

p_male_ribbons <- make_panel_plot(
  male_obj,
  panel_title = "Male",
  show_legend = FALSE
)

p_female_ribbons <- make_panel_plot(
  female_obj,
  panel_title = "Female",
  show_legend = TRUE
)

p_sankey_two_panels <- (
  p_male_ribbons / p_female_ribbons
) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    text = element_text(family = PATHWAY_FONT)
  )

p_sankey_two_panels

# Figure 5 uses the same Arial PNG/PDF devices as all other figures.
save_pathway_plot <- function(p, stem, height_mm = 225, width_mm = 180, dpi = 600) {
  save_main_gg(p, stem, height_mm, width_mm, dpi)
}

save_pathway_plot(
  p_sankey_two_panels,
  file.path(fig_main_dir, "Pathway_2_pans"),
  height_mm = 225
)

# 10. Figure 1 and supplementary preprocessing figures (graphics only).
export_preprocessing_figures <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  # Always use the exact saved analysis matrix for the After RUV panel.
  after <- as.matrix(readRDS(file.path(data_dir, "proteomic_cleaned_expr.rds")))
  meta <- readRDS(file.path(data_dir, "proteomic_nako_metadata_aligned.rds"))
  meta <- meta[match(colnames(after), meta$Sample_ID), , drop = FALSE]
  stopifnot(identical(as.character(meta$Sample_ID), colnames(after)))
  sex_plot <- factor(meta$sex, levels = c(1, 2), labels = c("male", "female"))
  age_plot <- factor(meta$ageclass, levels = 1:5,
                     labels = c("20-29", "30-39", "40-49", "50-59", "60-69"))
  
  before <- get0("uq_norm_log2", envir = parent.env(environment()),
                 inherits = TRUE, ifnotfound = NULL)
  if (!is.null(before)) {
    before <- as.matrix(before)
    if (!all(rownames(after) %in% rownames(before)) ||
        !all(colnames(after) %in% colnames(before))) {
      stop("The in-memory pre-RUV matrix does not match the saved analysis data.")
    }
    before <- before[rownames(after), colnames(after), drop = FALSE]
  } else {
    message("Reconstructing the pre-RUV matrix for Figure 1 only; analysis data are unchanged.")
    if (!requireNamespace("EDASeq", quietly = TRUE)) {
      stop("EDASeq is needed to reconstruct the pre-RUV plotting matrix.")
    }
    pre_path <- file.path(data_dir, "proteomic_preMICE_expr.rds")
    if (!file.exists(pre_path)) {
      stop("Run 01_preprocessing.R first: the preMICE checkpoint is missing.")
    }
    pre <- as.matrix(readRDS(pre_path))
    pre_meta <- readRDS(file.path(data_dir, "proteomic_preMICE_metadata.rds"))
    pre_meta <- pre_meta[match(colnames(pre), pre_meta$Sample_ID), , drop = FALSE]
    stopifnot(identical(as.character(pre_meta$Sample_ID), colnames(pre)))
    sex_before <- factor(pre_meta$sex, levels = c(1, 2), labels = c("male", "female"))
    sxp_plot <- t(pre)
    idx_by_sex_plot <- split(seq_len(nrow(sxp_plot)), sex_before)
    imp_plot <- lapply(names(idx_by_sex_plot), function(g) {
      idx <- idx_by_sex_plot[[g]]
      blk <- as.data.frame(sxp_plot[idx, , drop = FALSE])
      set.seed(123)
      mf_args <- list(
        xmis = blk, maxiter = 10, ntree = 500,
        mtry = floor(sqrt(ncol(blk))), replace = TRUE,
        nodesize = c(1, 5), parallelize = "no", verbose = TRUE
      )
      if ("backend" %in% names(formals(missForest::missForest))) {
        mf_args$backend <- "randomForest"
      }
      do.call(missForest::missForest, mf_args)$ximp
    })
    names(imp_plot) <- names(idx_by_sex_plot)
    imputed <- do.call(rbind, imp_plot)
    rownames(imputed) <- sub("^(male|female)\\.", "", rownames(imputed))
    imputed <- imputed[rownames(sxp_plot), , drop = FALSE]
    before <- log2(EDASeq::betweenLaneNormalization(t(as.matrix(imputed)),
                                                    which = "upper") + 1)
    stopifnot(all(rownames(after) %in% rownames(before)),
              all(colnames(after) %in% colnames(before)))
    before <- before[rownames(after), colnames(after), drop = FALSE]
  }
  
  # Identical standardisation/PCA to plot_pca() in 01_preprocessing.R.
  pca_panel <- function(M, main) {
    keep <- apply(M, 1, stats::var) > 0
    X <- t(scale(t(M[keep, , drop = FALSE]), center = TRUE, scale = TRUE))
    p <- prcomp(t(X))
    ve <- 100 * p$sdev^2 / sum(p$sdev^2)
    cols <- ifelse(sex_plot == "male", "#009E73", "#E69F00")
    plot(p$x[, 1], p$x[, 2], pch = 20 + as.integer(age_plot),
         col = "black", bg = cols, cex = 1.2, lwd = 1,
         xlab = sprintf("PC1 (%.1f%%)", ve[1]),
         ylab = sprintf("PC2 (%.1f%%)", ve[2]), main = main,
         cex.main = 1.25, cex.lab = 1.05, cex.axis = 0.95)
  }
  draw_pca <- function() {
    op <- par(no.readonly = TRUE)
    on.exit({ layout(1); par(op) }, add = TRUE)
    layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
           heights = c(4, 1.1))
    par(mar = c(4.6, 4.6, 3.4, 1.2), family = PLOT_FONT, cex = 1)
    pca_panel(before, "Before RUV")
    mtext("A", side = 3, line = 1.2, adj = 0, font = 2, cex = 1.5)
    pca_panel(after, "After RUV")
    mtext("B", side = 3, line = 1.2, adj = 0, font = 2, cex = 1.5)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    legend(x = 0.14, y = 0.45, xjust = 0.5, yjust = 0.5,
           legend = c("male", "female"), pt.bg = c("#009E73", "#E69F00"),
           pch = 21, col = "black", pt.cex = 1.2, cex = 0.95,
           horiz = TRUE, bty = "n", xpd = NA)
    # A separate title provides a deliberate blank gap above the symbols.
    text(0.64, 0.66, "Age category (years)", cex = 1.05)
    legend(x = 0.64, y = 0.45, xjust = 0.5, yjust = 0.5,
           legend = levels(age_plot), pch = 21:(20 + nlevels(age_plot)),
           col = "black", pt.bg = "grey85", pt.cex = 1.2, cex = 0.95,
           horiz = TRUE, bty = "n", xpd = NA)
  }
  save_main_base(draw_pca, file.path(fig_main_dir, "plot_pca_ruv"),
                 height_mm = 110)
  
  # The outlier replacement does not alter the NA mask. Thus the original
  # missingness curve can be redrawn directly from the raw CSV without rerunning QC.
  raw_plot <- read.csv(file.path(data_dir, "Nako_pdata.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  raw_values <- as.matrix(data.frame(lapply(raw_plot[, -1, drop = FALSE],
                                            function(x) as.numeric(gsub(",", ".", x, fixed = FALSE))), check.names = FALSE))
  thr_plot <- seq(0, 0.60, by = 0.005)
  miss_plot <- rowMeans(is.na(raw_values))
  retained_plot <- vapply(thr_plot, function(t) sum(miss_plot <= t), numeric(1))
  xn <- (thr_plot - min(thr_plot)) / diff(range(thr_plot))
  yn <- (retained_plot - min(retained_plot)) / diff(range(retained_plot))
  elbow_idx <- which.max(yn - xn)
  draw_elbow_tnr <- function() {
    par(mar = c(4, 4, 2, 1), family = PLOT_FONT)
    plot(thr_plot * 100, retained_plot, type = "l",
         xlab = "% NA threshold", ylab = "Proteins retained",
         main = "Global missingness (elbow)")
    abline(v = 100 * thr_plot[elbow_idx], lty = 2)
    points(100 * thr_plot[elbow_idx], retained_plot[elbow_idx], pch = 19)
    text(100 * thr_plot[elbow_idx], retained_plot[elbow_idx],
         labels = sprintf(" elbow ≈ %.1f%%", 100 * thr_plot[elbow_idx]),
         pos = 4, cex = 0.9)
  }
  save_base_plot(draw_elbow_tnr, file.path(fig_supp_dir, "elbow_plot"),
                 width_px = 1565, height_px = 862)
  invisible(NULL)
}

# Complete Figure S1: text transcribed from flowchart_preprocessing_workflow_complete.png.
# This is a fixed diagram of the supplied workflow, not another analysis.
# It uses grid (included with R); no flowchart package is needed.
draw_complete_preprocessing_flowchart <- function() {
  headings <- c(
    "Input data", "Sample-level QC", "Outlier handling", "Metadata alignment",
    "Missingness assessment", "Protein filtering", "Missing value imputation",
    "Between-sample normalization", "Log2 transform for analysis",
    "Remove unwanted technical variation", "Post-normalization QC", "Output"
  )
  details <- c(
    "Raw plasma protein intensities (488 proteins) from selected NAKO subsample (117 participants)",
    "Per-sample missingness check + PCA outlier screening (|z| > 4 on PC1/PC2; no samples excluded)",
    "Per protein: censor values outside mean ± 2.5 SD and impute using a censored log-normal model\n(1141 value replaced = 2.342% of total values)",
    "Align samples to NAKO metadata (sex; 5 age groups; bmi, smoke, education, alcohol)",
    "Test associations between missingness and key covariates (sex, age group);\nno significant associations detected → consistent with MCAR",
    "Filter proteins by missingness using a data-driven elbow threshold (Kneedle); retained 370 proteins",
    "Impute remaining missing values (<1% after filtering) using missForest within sex strata",
    "Upper-quartile normalization",
    "log2(normalized + 1)",
    "RUV adjustment to remove technical variation while preserving age and sex signal (negative controls:\nproteins minimally associated with age in both sexes; k tuned to preserve sex separation)",
    "RLE + PCA diagnostics (improved preservation of biological sex separation)",
    "Cleaned expression matrix (370 proteins × 117 participants; 60 males, 57 females) + aligned metadata"
  )
  grid::grid.newpage()
  grid::grid.text("Preprocessing Workflow Flow Chart", x = 0.5, y = 0.965,
                  gp = grid::gpar(fontfamily = PLOT_FONT, fontsize = 18, fontface = "bold"))
  # Measure the text on the active device so that both PNG and PDF get safe boxes.
  body_gp <- grid::gpar(fontfamily = PLOT_FONT, fontsize = 11.5, lineheight = 1.15)
  head_gp <- grid::gpar(fontfamily = PLOT_FONT, fontsize = 12, fontface = "bold")
  body_grobs <- lapply(details, function(x) grid::textGrob(x, gp = body_gp))
  title_grobs <- lapply(headings, function(x) grid::textGrob(x, gp = head_gp))
  body_h <- vapply(body_grobs, function(g) grid::convertHeight(grid::grobHeight(g), "npc", TRUE), numeric(1))
  head_h <- vapply(title_grobs, function(g) grid::convertHeight(grid::grobHeight(g), "npc", TRUE), numeric(1))
  text_w <- max(vapply(c(body_grobs, title_grobs), function(g)
    grid::convertWidth(grid::grobWidth(g), "npc", TRUE), numeric(1)))
  box_w <- max(0.90, text_w + 0.025)
  if (box_w > 0.97) stop("Flowchart text is too wide: increase its export width.")
  heights <- body_h + head_h + 0.016
  gap <- (0.88 - sum(heights)) / (length(headings) - 1)
  if (gap < 0.01) stop("Flowchart text is too tall: increase its export height.")
  top <- 0.925
  for (i in seq_along(headings)) {
    bottom <- top - heights[i]
    grid::grid.rect(x = 0.5, y = (top + bottom) / 2, width = box_w, height = heights[i],
                    gp = grid::gpar(fill = "white", col = "black", lwd = 1.2))
    grid::grid.text(headings[i], x = 0.5, y = top - 0.005 - head_h[i] / 2, gp = head_gp)
    grid::grid.text(details[i], x = 0.5, y = bottom + 0.005 + body_h[i] / 2, gp = body_gp)
    if (i < length(headings)) {
      grid::grid.lines(x = c(0.5, 0.5),
                       y = c(bottom - 0.002, bottom - gap + 0.002),
                       arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "open"),
                       gp = grid::gpar(col = "black", lwd = 0.8))
    }
    top <- bottom - gap
  }
}

export_preprocessing_figures()
save_base_plot(draw_complete_preprocessing_flowchart,
               file.path(fig_supp_dir, "flowchart_preprocessing_workflow_complete"),
               width_px = 1800, height_px = 2300, dpi = 180)
# Keep the existing standard Figure S1 filename too, using the same complete diagram.
for (.ext in c("png", "pdf")) {
  .from <- file.path(fig_supp_dir, paste0("flowchart_preprocessing_workflow_complete.", .ext))
  .to <- file.path(fig_supp_dir, paste0("flowchart_preprocessing_workflow.", .ext))
  if (!file.copy(.from, .to, overwrite = TRUE)) stop("Could not save Figure S1: ", .to)
}

# Compact run summary.
summary_file <- file.path(results_dir, "main_analysis_summary.txt")
sink(summary_file)
cat("Final sample size:", ncol(expr_clean), "\n")
cat("Final protein count:", nrow(expr_clean), "\n")
cat("Stable proteins - male:", length(stable_ids_m), "\n")
cat("Stable proteins - female:", length(stable_ids_f), "\n")
cat("AgeProt R2 male:", sprintf("%.3f", summary(lm(AgeProt_score ~ age10, data = res_male$data))$r.squared), "\n")
cat("AgeProt R2 female:", sprintf("%.3f", summary(lm(AgeProt_score ~ age10, data = res_female$data))$r.squared), "\n")
sink()

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_analysis.txt"))
message("Main analysis complete.")
