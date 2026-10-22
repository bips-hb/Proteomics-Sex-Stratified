# 01_preprocessing.R


suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(data.table)
  library(matrixStats)
  library(lnormimp)
  library(EDASeq)
  library(missForest)
  library(PhosR)
  library(limma)
  library(cluster)
  library(MatrixGenerics)
  library(ggplot2)
})

here::i_am("analysis/code/01_preprocessing.R")

set.seed(123)

cat("Reproducibility-critical package versions:\n")
cat("  R:", R.version.string, "\n")
cat("  missForest:", as.character(utils::packageVersion("missForest")), "\n")
cat("  randomForest:", as.character(utils::packageVersion("randomForest")), "\n")
cat("  PhosR:", as.character(utils::packageVersion("PhosR")), "\n")
cat("  EDASeq:", as.character(utils::packageVersion("EDASeq")), "\n")

data_dir <- here("analysis", "data")
fig_main_dir <- here("analysis", "results", "figures_main")
fig_supp_dir <- here("analysis", "results", "figures_supplementary")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp_dir, recursive = TRUE, showWarnings = FALSE)

save_base_plot <- function(draw_fun, stem, width_px, height_px, dpi = 96) {
  png(file.path(dirname(stem), paste0(basename(stem), ".png")),
      width = width_px, height = height_px, res = dpi)
  draw_fun()
  dev.off()
  pdf(file.path(dirname(stem), paste0(basename(stem), ".pdf")),
      width = width_px / dpi, height = height_px / dpi, useDingbats = FALSE)
  draw_fun()
  dev.off()
}

# 1. Load raw proteomics
f <- here("analysis", "data", "Nako_pdata.csv")
proteomic_data <- fread(f, data.table = FALSE)

cat(sprintf("Number of proteins: %d\nNumber of subjects: %d\n",
            dim(proteomic_data)[1], dim(proteomic_data)[2]-1))

prot_ids <- proteomic_data$proteins
X_raw_chr <- proteomic_data[, -1, drop = FALSE]
to_num <- function(x) as.numeric(gsub(",", ".", x, fixed = FALSE))
X_raw_num <- as.data.frame(lapply(X_raw_chr, to_num), check.names = FALSE)
rownames(X_raw_num) <- prot_ids

expr <- as.matrix(X_raw_num)
rownames(expr) <- proteomic_data$proteins
log_expr <- log2(expr + 1)

log_expr_pca <- log_expr

keep_rows <- rowSums(!is.na(log_expr_pca)) >= 2
log_expr_pca <- log_expr_pca[keep_rows, , drop = FALSE]

row_med <- matrixStats::rowMedians(log_expr_pca, na.rm = TRUE)
na_idx <- which(is.na(log_expr_pca), arr.ind = TRUE)
if (nrow(na_idx) > 0) log_expr_pca[na_idx] <- row_med[na_idx[,1]]

non_zero_var <- matrixStats::rowVars(log_expr_pca) > 0
log_expr_pca <- log_expr_pca[non_zero_var, , drop = FALSE]

pca <- prcomp(t(log_expr_pca), center = TRUE, scale. = TRUE)

z_pc1 <- as.numeric(scale(pca$x[,1]))
z_pc2 <- as.numeric(scale(pca$x[,2]))

thr_z <- 4
flag_pca_simple <- (abs(z_pc1) > thr_z) | (abs(z_pc2) > thr_z)
suspects_pca_simple <- colnames(log_expr_pca)[flag_pca_simple]

if (length(suspects_pca_simple)) {
  message("PCA outliers (|z| > ", thr_z, " on PC1 or PC2): ",
          paste(suspects_pca_simple, collapse = ", "))
} else {
  message("No PCA outliers at |z| > ", thr_z, ".")
}

# 2. Outlier handling
set.seed(123)

# Derives LC/UC using the mean ± z·SD rule
derive_cutoffs <- function(x, z = 2.5) {
  ok <- is.finite(x) & x > 0
  if (sum(ok) < 3L) return(c(NA_real_, NA_real_))
  m  <- mean(x[ok])
  s  <- stats::sd(x[ok])
  if (!is.finite(s) || s <= 0) return(c(NA_real_, NA_real_))
  lc <- max(m - z * s, .Machine$double.xmin)
  uc <- m + z * s
  c(lc, uc)
}

# Row-wise Herbers-style imputation using lnormimp::lnormimp
impute_herbers_row <- function(x, z = 2.5, lc = NA_real_, uc = NA_real_) {
  ok <- is.finite(x) & x > 0
  if (sum(ok) < 3L) return(x)
  
  # If LC/UC are not provided, estimate via mean ± z·SD
  if (!is.finite(lc) || !is.finite(uc)) {
    cu <- derive_cutoffs(x[ok], z = z)
    lc <- cu[1]; uc <- cu[2]
  }
  if (!is.finite(lc) || !is.finite(uc) || lc >= uc) return(x)
  
  # Identify censored values outside [lc, uc]
  below  <- ok & (x < lc)
  above  <- ok & (x > uc)
  inside <- ok & !below & !above
  
  n_below <- sum(below)
  n_above <- sum(above)
  if ((n_below + n_above) == 0L) return(x)
  
  data_cens <- x[inside]
  
  # Fit censored log-normal model (MLE) and generate imputations (Herbers et al.)
  imp <- tryCatch(
    lnormimp::lnormimp(
      data_censored = data_cens,
      censn         = c(n_below, n_above),
      cutoff        = c(lc, uc)
    ),
    error = function(e) NULL
  )
  if (is.null(imp)) return(x)
  
  # Reconstruct row: the last (n_below + n_above) entries are the imputed draws
  repl <- utils::tail(imp, n_below + n_above)
  x[below] <- repl[seq_len(n_below)]
  x[above] <- repl[n_below + seq_len(n_above)]
  x
}

# If LLOQ/ULOQ are supplied, they are used; otherwise LC/UC are estimated via mean ± z·SD.
impute_herbers <- function(expr, z = 2.5, LLOQ = NULL, ULOQ = NULL) {
  stopifnot(is.matrix(expr))
  has_limits <- !is.null(LLOQ) && !is.null(ULOQ) &&
    length(LLOQ) == nrow(expr) && length(ULOQ) == nrow(expr)
  out <- expr
  for (i in seq_len(nrow(expr))) {
    lc <- if (has_limits) LLOQ[i] else NA_real_
    uc <- if (has_limits) ULOQ[i] else NA_real_
    out[i, ] <- impute_herbers_row(expr[i, ], z = z, lc = lc, uc = uc)
  }
  dimnames(out) <- dimnames(expr)
  out
}

# --- Run imputation (z = 2.5 SD) ---
expr_outlier_fixed <- impute_herbers(expr, z = 2.5)

# --- Quick stats: count imputed values ---
count_replaced <- function(x, z = 2.5, lc = NA_real_, uc = NA_real_) {
  ok <- is.finite(x) & x > 0
  if (sum(ok) < 3L) return(0L)
  if (!is.finite(lc) || !is.finite(uc)) {
    cu <- derive_cutoffs(x[ok], z = z)
    if (any(!is.finite(cu))) return(0L)
    lc <- cu[1]; uc <- cu[2]
  }
  sum(ok & (x < lc | x > uc))
}

n_replaced <- sum(apply(expr, 1, count_replaced, z = 2.5))
n_total    <- sum(is.finite(expr))
message(sprintf(
  "Herbers-style censored log-normal (z=%.1f SD): %d replaced (%.3f%% of finite).",
  2.5, n_replaced, 100 * n_replaced / n_total
))

# 3. Align NAKO metadata and assess missingness
nako <- readRDS(here("analysis", "data", "nako_base.rds"))   # contains sex, age, ageclass, group, Sample_ID

sample_ids <- colnames(expr_outlier_fixed)
meta_aligned <- nako %>%
  dplyr::filter(Sample_ID %in% sample_ids) %>%
  dplyr::arrange(match(Sample_ID, sample_ids))

stopifnot(identical(meta_aligned$Sample_ID, sample_ids))

# Extract covariates
sex_vec      <- factor(meta_aligned$sex, levels = c(1,2), labels = c("male","female"))
age_num      <- meta_aligned$age
ageclass_fac <- factor(meta_aligned$ageclass, levels = 1:5,
                       labels = c("20-<30","30-<40","40-<50","50-<60",">=60"))
grps         <- factor(meta_aligned$group)  # e.g., "m:20-<30", ... "f:>=60"

X    <- expr_outlier_fixed                    # matrix: proteins x samples
M    <- is.na(X)                              # missingness matrix (TRUE = NA)
stopifnot(identical(colnames(X), meta_aligned$Sample_ID))

# --- Sample-level tests (NA% per sample vs covariates)
na_perc_samp <- 100 * colMeans(M)
# Sex: Wilcoxon
p_wilcox_sex <- tryCatch(wilcox.test(na_perc_samp ~ sex_vec)$p.value, error = function(e) NA_real_)
# Age: Spearman
p_spear_age  <- suppressWarnings(cor.test(na_perc_samp, age_num, method = "spearman")$p.value)
# Age class: Kruskal–Wallis
p_kw_agecl   <- tryCatch(kruskal.test(na_perc_samp ~ ageclass_fac)$p.value, error = function(e) NA_real_)
# Group (sex×ageclass): Kruskal–Wallis
p_kw_group   <- tryCatch(kruskal.test(na_perc_samp ~ grps)$p.value, error = function(e) NA_real_)

cat(sprintf(
  "Sample-level NA%% tests: Sex (Wilcoxon) p=%.3g | Age (Spearman) p=%.3g | \n AgeClass (KW) p=%.3g | Group (KW) p=%.3g\n",
  p_wilcox_sex, p_spear_age, p_kw_agecl, p_kw_group))

# --- Protein-level: logistic GLMs for missingness
X <- expr_outlier_fixed
M <- is.na(X)                          # 1 = missing
prot <- rownames(X)

fit_one <- function(y) {
  # y: 0/1 vector (length = n samples)
  if (length(unique(y)) < 2) return(rep(NA_real_, 3))  # all 0 or all 1
  # full vs null (global LRT)
  fit0 <- glm(y ~ 1, family = binomial())
  fit1 <- glm(y ~ sex_vec + ageclass_fac, family = binomial())
  p_global <- tryCatch(anova(fit0, fit1, test = "Chisq")$`Pr(>Chi)`[2], error = function(e) NA_real_)
  dd <- tryCatch(drop1(fit1, test = "Chisq"), error = function(e) NULL)
  if (is.null(dd)) return(c(p_global, NA, NA))
  p_sex      <- dd["sex_vec",      "Pr(>Chi)"]
  p_ageclass <- dd["ageclass_fac", "Pr(>Chi)"]
  c(p_global, p_sex, p_ageclass)
}

P <- t(apply(M, 1, fit_one))
colnames(P) <- c("p_global","p_sex","p_ageclass")

res_prot <- data.frame(protein = prot, P, row.names = NULL, check.names = FALSE)

# BH correction (per colonna)
for (cc in c("p_global","p_sex","p_ageclass")) {
  res_prot[[paste0(cc,"_BH")]] <- p.adjust(res_prot[[cc]], method = "BH")
}

res_prot <- res_prot[order(res_prot$p_global), ]

# Quick summary + top hits
cat("Significant at FDR<5% — global:", sum(res_prot$p_global_BH < 0.05, na.rm = TRUE),
    "sex:",      sum(res_prot$p_sex_BH < 0.05, na.rm = TRUE),
    "ageclass:", sum(res_prot$p_ageclass_BH < 0.05, na.rm = TRUE), "\n")

head(res_prot, 10)

# ============================================================
# Protein-level: logistic GLMs for missingness vs
# BMI / education / smoke / alcohol  
# ============================================================

X    <- expr_outlier_fixed
M    <- is.na(as.matrix(X))            # TRUE = missing
prot <- rownames(X)

stopifnot(identical(colnames(X), meta_aligned$Sample_ID))

# --------- Helper: find a column name robustly by regex patterns ----------
find_col <- function(df, patterns) {
  nm <- names(df)
  hit <- unique(unlist(lapply(patterns, function(p) grep(p, nm, ignore.case = TRUE, value = TRUE))))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# --------- Identify covariate columns (edit patterns if needed) ----------
bmi_col     <- find_col(meta_aligned, c("^bmi$", "body.?mass"))
educ_col    <- find_col(meta_aligned, c("^education$", "^educ", "school", "bildung"))
smoke_col   <- find_col(meta_aligned, c("^smoke$", "smok", "rauch"))
alcohol_col <- find_col(meta_aligned, c("^alcohol$", "alc", "alkohol"))

cat("Detected columns:\n",
    "  BMI     ->", bmi_col, "\n",
    "  Educ    ->", educ_col, "\n",
    "  Smoke   ->", smoke_col, "\n",
    "  Alcohol ->", alcohol_col, "\n")

# Fail fast with a clear message if something wasn't found
if (any(is.na(c(bmi_col, educ_col, smoke_col, alcohol_col)))) {
  stop("Could not find one or more covariate columns in meta_aligned. ",
       "Run names(meta_aligned) and adjust the regex patterns in find_col().")
}

# --------- Extract / coerce covariates ----------
bmi_vec   <- suppressWarnings(as.numeric(meta_aligned[[bmi_col]]))
educ_fac  <- as.factor(meta_aligned[[educ_col]])
smoke_fac <- as.factor(meta_aligned[[smoke_col]])

# Alcohol: try numeric; if too sparse, treat as categorical
alcohol_try <- suppressWarnings(as.numeric(meta_aligned[[alcohol_col]]))
if (sum(!is.na(alcohol_try)) < 0.5 * length(alcohol_try)) {
  alcohol_var <- as.factor(meta_aligned[[alcohol_col]])
  alcohol_is_factor <- TRUE
} else {
  alcohol_var <- alcohol_try
  alcohol_is_factor <- FALSE
}
cat("Alcohol treated as:", ifelse(alcohol_is_factor, "factor", "numeric"), "\n")

# --------- Per-protein GLM function ----------
fit_one_covs <- function(y) {
  out <- rep(NA_real_, 5)
  names(out) <- c("p_global","p_bmi","p_educ","p_smoke","p_alcohol")
  
  y <- as.integer(y)
  if (length(unique(y)) < 2) return(unname(out))  # all observed or all missing
  
  dat <- data.frame(
    y = y,
    bmi = bmi_vec,
    educ = educ_fac,
    smoke = smoke_fac,
    alcohol = alcohol_var
  )
  
  # Complete-case on covariates (so we test missingness vs covariates cleanly)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < 10) return(unname(out))
  if (length(unique(dat$y)) < 2) return(unname(out))
  
  # Guard: factors must have >=2 levels after complete-case
  if (is.factor(dat$educ)   && nlevels(droplevels(dat$educ))   < 2) return(unname(out))
  if (is.factor(dat$smoke)  && nlevels(droplevels(dat$smoke))  < 2) return(unname(out))
  if (is.factor(dat$alcohol) && nlevels(droplevels(dat$alcohol)) < 2) return(unname(out))
  
  fit0 <- tryCatch(glm(y ~ 1, data = dat, family = binomial()),
                   error = function(e) NULL)
  fit1 <- tryCatch(glm(y ~ bmi + educ + smoke + alcohol, data = dat, family = binomial()),
                   error = function(e) NULL)
  
  if (is.null(fit0) || is.null(fit1)) return(unname(out))
  
  # Global likelihood-ratio test (null vs full)
  out["p_global"] <- tryCatch(anova(fit0, fit1, test = "Chisq")$`Pr(>Chi)`[2],
                              error = function(e) NA_real_)
  
  # BMI: Wald p-value (single coefficient) if present
  co <- tryCatch(summary(fit1)$coefficients, error = function(e) NULL)
  if (!is.null(co) && "bmi" %in% rownames(co)) out["p_bmi"] <- co["bmi","Pr(>|z|)"]
  
  # For multi-df terms (factors), use term-wise LRTs
  out["p_educ"] <- tryCatch({
    fit_educ0 <- update(fit1, . ~ . - educ)
    anova(fit_educ0, fit1, test = "Chisq")$`Pr(>Chi)`[2]
  }, error = function(e) NA_real_)
  
  out["p_smoke"] <- tryCatch({
    fit_smoke0 <- update(fit1, . ~ . - smoke)
    anova(fit_smoke0, fit1, test = "Chisq")$`Pr(>Chi)`[2]
  }, error = function(e) NA_real_)
  
  out["p_alcohol"] <- tryCatch({
    fit_alc0 <- update(fit1, . ~ . - alcohol)
    anova(fit_alc0, fit1, test = "Chisq")$`Pr(>Chi)`[2]
  }, error = function(e) NA_real_)
  
  unname(out)
}

# --------- Run per-protein ----------
P2 <- t(vapply(seq_len(nrow(M)),
               function(i) fit_one_covs(M[i, ]),
               FUN.VALUE = numeric(5)))

colnames(P2) <- c("p_global","p_bmi","p_educ","p_smoke","p_alcohol")

res_prot_covs <- data.frame(protein = prot, P2, row.names = NULL, check.names = FALSE)

# --------- Multiple testing correction (BH/FDR) per p-value column ----------
for (cc in colnames(P2)) {
  res_prot_covs[[paste0(cc, "_BH")]] <- p.adjust(res_prot_covs[[cc]], method = "BH")
}

# Sort by global p-value and print summary
res_prot_covs <- res_prot_covs[order(res_prot_covs$p_global), ]

cat("Significant at FDR<5% — global:",
    sum(res_prot_covs$p_global_BH < 0.05, na.rm = TRUE),
    "bmi:",     sum(res_prot_covs$p_bmi_BH < 0.05, na.rm = TRUE),
    "educ:",    sum(res_prot_covs$p_educ_BH < 0.05, na.rm = TRUE),
    "smoke:",   sum(res_prot_covs$p_smoke_BH < 0.05, na.rm = TRUE),
    "alcohol:", sum(res_prot_covs$p_alcohol_BH < 0.05, na.rm = TRUE),
    "\n")

head(res_prot_covs, 10)
# ============================================================

# 4. Protein missingness filtering
# --- NA filtering: global-only elbow (Kneedle) ---

X <- expr_outlier_fixed                       # proteins x samples (raw scale after outlier fix)
thr_seq <- seq(0, 0.60, by = 0.005)          # candidate NA thresholds (proportion)

global_miss     <- rowMeans(is.na(X))
n_keep_global   <- sapply(thr_seq, function(t) sum(global_miss <= t))

# Kneedle-style elbow: max vertical gap from min–max scaled diagonal
elbow_kneedle <- function(x, y) {
  xn <- (x - min(x)) / diff(range(x))
  yn <- (y - min(y)) / diff(range(y))
  i  <- which.max(yn - xn)
  list(idx = i, thr = x[i], kept = y[i])
}
e_glob <- elbow_kneedle(thr_seq, n_keep_global)

draw_elbow <- function() {
  par(mar = c(4,4,2,1))
  plot(thr_seq * 100, n_keep_global, type = "l",
       xlab = "% NA threshold", ylab = "Proteins retained",
       main = "Global missingness (elbow)")
  abline(v = 100 * e_glob$thr, lty = 2)
  points(100 * e_glob$thr, e_glob$kept, pch = 19)
  text(100 * e_glob$thr, e_glob$kept,
       labels = sprintf(" elbow ≈ %.1f%%", 100 * e_glob$thr),
       pos = 4, cex = 0.9)
}
# Apply global threshold
thr_global  <- e_glob$thr
keep_global <- global_miss <= thr_global
expr_filtered <- X[keep_global, , drop = FALSE]

# Summary
n <- nrow(X)
rm_global <- sum(!keep_global)
pct <- function(k) 100 * k / n
cat(sprintf("Global NA filtering: %.1f%%. Removed %d (%.2f%%), Kept %d (%.2f%%)\n",
            100 * thr_global, rm_global, pct(rm_global),
            sum(keep_global), pct(sum(keep_global))))

save_base_plot(
  draw_elbow,
  file.path(fig_supp_dir, "elbow_plot"),
  width_px = 1565, height_px = 862
)

# Checkpoint immediately before missing-value imputation.
saveRDS(expr_filtered, file.path(data_dir, "proteomic_preMICE_expr.rds"))
saveRDS(meta_aligned, file.path(data_dir, "proteomic_preMICE_metadata.rds"))

# 5. MissForest imputation within sex
# 1) Work on samples x proteins
sxp <- t(expr_filtered)                  # rows = samples, cols = proteins
stopifnot(nrow(sxp) == length(sex_vec))  # alignment check

n_na_before <- sum(is.na(sxp))
tot_cells   <- length(sxp)

# 2) Split by sex, impute each stratum
idx_by_sex <- split(seq_len(nrow(sxp)), sex_vec)

imp_list <- lapply(names(idx_by_sex), function(g) {
  
  idx <- idx_by_sex[[g]]
  blk <- as.data.frame(sxp[idx, , drop = FALSE])
  
  set.seed(123)
  
  # Reproduce missForest 1.5 + randomForest semantics explicitly.
  # In missForest 1.5 the effective default was nodesize = 1 for numeric
  # responses and 5 for factors. missForest >= 1.6 changed that default.
  mf_args <- list(
    xmis = blk,
    maxiter = 10,
    ntree = 500,
    mtry = floor(sqrt(ncol(blk))),
    replace = TRUE,
    nodesize = c(1, 5),
    parallelize = "no",
    verbose = TRUE
  )
  
  # missForest 1.5 has no 'backend' argument; >=1.6 does.
  # When present, force the historical randomForest implementation.
  if ("backend" %in% names(formals(missForest::missForest))) {
    mf_args$backend <- "randomForest"
  }
  
  mf <- do.call(missForest::missForest, mf_args)
  mf$ximp
})
names(imp_list) <- names(idx_by_sex)

# 3) Rebind strata back in the original sample order
sxp_imp <- do.call(rbind, imp_list)
# safer: only remove "male." / "female." prefix from row names
rownames(sxp_imp) <- sub("^(male|female)\\.", "", rownames(sxp_imp))
sxp_imp <- sxp_imp[rownames(sxp), , drop = FALSE]  # original sample order
stopifnot(identical(rownames(sxp_imp), rownames(sxp)))

# 4) Back to proteins x samples, restore dimnames
expr_imputed <- t(as.matrix(sxp_imp))
rownames(expr_imputed) <- rownames(expr_filtered)
colnames(expr_imputed) <- colnames(expr_filtered)

# Quick sanity check + percent imputed
n_na_after <- sum(is.na(expr_imputed))
n_imputed  <- n_na_before - n_na_after
perc_imputed <- 100 * n_imputed / tot_cells

cat("Missing values before imputation:", n_na_before, "\n")
cat("Remaining NAs after missForest (should be 0):", n_na_after, "\n")
cat(sprintf("Imputed values: %d (%.3f%% of all cells)\n", n_imputed, perc_imputed))

# 6. UQ normalization, log2 and RUV
# Apply between-lane (sample-wise) normalization using the upper quantile (on non-log scale)
uq_norm_raw <- EDASeq::betweenLaneNormalization(expr_imputed, which = "upper")

uq_norm_log2 <- log2(uq_norm_raw + 1)

## ----- Control selection (independent of k except for min_controls) -----
target_prop <- 0.10

# Exclude sentinel(s)
pzp_name <- grep("^(P20742|PZP(_HUMAN)?)$", rownames(uq_norm_log2),
                 ignore.case = TRUE, value = TRUE)

# Age-association per sex (to rank negative controls)
pvals_by_sex <- lapply(levels(sex_vec), function(sx) {
  idx <- which(sex_vec == sx)
  Y <- uq_norm_log2[, idx, drop = FALSE]
  age_s <- age_num[idx]
  fit <- eBayes(lmFit(Y, model.matrix(~ age_s)))
  data.frame(protein = rownames(Y), pval = fit$p.value[, "age_s"])
})
names(pvals_by_sex) <- levels(sex_vec)

tab <- Reduce(function(a, b) merge(a, b, by = "protein", suffixes = c("_m","_f")), pvals_by_sex)
tab <- subset(tab, !(protein %in% pzp_name))
tab$score <- pmin(tab$pval_m, tab$pval_f, na.rm = TRUE) # bigger = less associated with age (for both sex)
tab <- tab[order(-tab$score), ]

n_target <- round(target_prop * nrow(uq_norm_log2))  # base target controls (fixed)

## ----- Helper: average silhouette on PZP -----
.mean_silhouette <- function(x, groups) {
  # x: numerical vector (es. PZP expression over all the samples)
  # groups: group factor (es. sex_vec)
  g <- droplevels(base::as.factor(groups))
  if (nlevels(g) < 2 || length(x) < 3 || any(table(g) < 2)) return(NA_real_)
  d <- dist(matrix(x, ncol = 1))  # distanza euclidea 1D
  sil <- suppressWarnings(cluster::silhouette(as.integer(g), d))
  mean(sil[, "sil_width"], na.rm = TRUE)
}

## ----- Grid search over k using PZP sex silhouette as criterion -----
k_grid <- 0:10

run_ruv_k <- function(k) {
  min_controls <- max(100, k * 10) # has to be quite big
  n_pick <- min(max(n_target, min_controls), nrow(tab))
  controls <- head(tab$protein, n_pick)
  ctl_idx  <- which(rownames(uq_norm_log2) %in% controls)
  
  ruv_list <- lapply(levels(sex_vec), function(sx) {
    idx <- which(sex_vec == sx)
    Y   <- uq_norm_log2[, idx, drop = FALSE]
    M   <- model.matrix(~ 0 + factor(ageclass_fac[idx]))
    PhosR::RUVphospho(Y, M = M, k = k, ctl = ctl_idx)
  })
  out <- do.call(cbind, ruv_list)
  out[, colnames(uq_norm_log2), drop = FALSE]
}

scores <- sapply(k_grid, function(k) {
  U <- run_ruv_k(k)
  x <- as.numeric(U[pzp_name, ])
  .mean_silhouette(x, sex_vec)
})

best_idx <- which.max(scores)
k_grid_best <- k_grid[best_idx]

k_ruv <- 8L
n_pick <- 100L

cat(sprintf(
  "Sentinel grid-search best k = %d (silhouette = %.3f); manuscript k fixed at %d.\n",
  k_grid_best, scores[best_idx], k_ruv
))

if (k_grid_best != k_ruv) {
  warning(
    "Current package/environment selects k=", k_grid_best,
    " by the diagnostic grid search; manuscript results used k=8."
  )
}

if (nrow(tab) < n_pick) {
  stop("Fewer than 100 candidate negative controls are available.")
}

controls <- head(tab$protein, n_pick)
ctl_idx  <- which(rownames(uq_norm_log2) %in% controls)
stopifnot(length(ctl_idx) == 100L)

ruv_list <- lapply(levels(sex_vec), function(sx) {
  idx <- which(sex_vec == sx)
  Y   <- uq_norm_log2[, idx, drop = FALSE]
  M   <- model.matrix(~ 0 + factor(ageclass_fac[idx]))
  PhosR::RUVphospho(Y, M = M, k = k_ruv, ctl = ctl_idx)
})
ruv_mat     <- do.call(cbind, ruv_list)
uq_ruv_log2 <- ruv_mat[, colnames(uq_norm_log2), drop = FALSE]

# Quick logs (silhouette su PZP post-RUV)
pzp_sil <- .mean_silhouette(as.numeric(uq_ruv_log2[pzp_name, ]), sex_vec)
cat(sprintf("RUV completed with k = %d and %d controls | post-RUV PZP sex silhouette = %.3f\n",
            k_ruv, length(ctl_idx), pzp_sil))
cat("Matrix dims post-RUV:", paste(dim(uq_ruv_log2), collapse = " x "), "\n")

# Manuscript reproduction checks.
stopifnot(
  nrow(uq_ruv_log2) == 370L,
  ncol(uq_ruv_log2) == 117L,
  length(ctl_idx) == 100L,
  k_ruv == 8L
)

.pca_variance <- function(M) {
  keep <- MatrixGenerics::rowVars(M) > 0
  X <- t(scale(t(M[keep, , drop = FALSE]), center = TRUE, scale = TRUE))
  p <- prcomp(t(X))
  100 * (p$sdev^2 / sum(p$sdev^2))
}

ve_pre  <- .pca_variance(uq_norm_log2)
ve_post <- .pca_variance(uq_ruv_log2)

cat(sprintf(
  "PCA variance: pre-RUV PC1/PC2 = %.1f/%.1f%% | post-RUV PC1/PC2 = %.1f/%.1f%%\n",
  ve_pre[1], ve_pre[2], ve_post[1], ve_post[2]
))

cat("\n================ MANUSCRIPT DIAGNOSTICS ================\n")

cat("Samples:", ncol(uq_norm_log2), "\n")
cat("Proteins:", nrow(uq_norm_log2), "\n")

cat(
  "Missing before MissForest:", n_na_before,
  "\nMissing after MissForest:", n_na_after, "\n"
)

cat(
  "Pre-RUV range:",
  min(uq_norm_log2, na.rm = TRUE),
  "to",
  max(uq_norm_log2, na.rm = TRUE),
  "\n"
)

cat(
  "Grid-search best k:",
  k_grid_best,
  "\n"
)

cat(
  "Grid-search silhouette:",
  scores[which(k_grid == k_grid_best)],
  "\n"
)

cat(
  "Forced manuscript k:",
  k_ruv,
  "\n"
)

cat(
  "Number of RUV controls:",
  length(ctl_idx),
  "\n"
)

cat(
  "Post-RUV PZP silhouette:",
  pzp_sil,
  "\n"
)

cat(
  "Post-RUV range:",
  min(uq_ruv_log2, na.rm = TRUE),
  "to",
  max(uq_ruv_log2, na.rm = TRUE),
  "\n"
)

cat(
  "Post-RUV PCA PC1/PC2:",
  round(ve_post[1], 3),
  round(ve_post[2], 3),
  "\n"
)

cat("=========================================================\n")

cat(sprintf(
  "Historical manuscript PCA diagnostic: post-RUV PC1/PC2 = %.1f/%.1f%%\n",
  ve_post[1], ve_post[2]
))

cat(sprintf(
  "Ranges:\n  expr (raw):       %s to %s\n  uq_norm_log2:     %s to %s\n  uq_ruv_log2:      %s to %s\n",
  min(expr, na.rm = TRUE),        max(expr, na.rm = TRUE),
  min(uq_norm_log2, na.rm = TRUE), max(uq_norm_log2, na.rm = TRUE),
  min(uq_ruv_log2, na.rm = TRUE),  max(uq_ruv_log2, na.rm = TRUE)
))

plot_pca <- function(M, main, sex_vec, ageclass_fac,
                     cex_pts   = 1.7,
                     lwd_pts   = 1.8,
                     cex_main  = 1.6,   
                     cex_lab   = 1.35,  
                     cex_axis  = 1.15,  
                     cex_A_B   = 1.7    
) {
  keep <- MatrixGenerics::rowVars(M) > 0
  X <- t(scale(t(M[keep, , drop = FALSE]), center = TRUE, scale = TRUE))
  p <- prcomp(t(X))
  ve <- 100 * (p$sdev^2 / sum(p$sdev^2))
  
  # Color-blind safe (Okabe–Ito)
  col_m <- "#009E73"  # blue-green
  col_f <- "#E69F00"  # orange
  cols  <- ifelse(sex_vec == "male", col_m, col_f)
  
  # Age class shapes: 21–25 filled symbols
  pchs <- 20 + as.integer(ageclass_fac)  # 21..25
  
  op_loc <- par(cex.main = cex_main, cex.lab = cex_lab, cex.axis = cex_axis)
  on.exit(par(op_loc), add = TRUE)
  
  plot(
    p$x[, 1], p$x[, 2],
    pch  = pchs,
    col  = "black",
    bg   = cols,
    cex  = cex_pts,
    lwd  = lwd_pts,
    xlab = sprintf("PC1 (%.1f%%)", ve[1]),
    ylab = sprintf("PC2 (%.1f%%)", ve[2]),
    main = main
  )
  
  invisible(list(col_m = col_m, col_f = col_f, cex_A_B = cex_A_B))
}

# 7. Manuscript figures and processed outputs
draw_pca_ruv <- function() {
  # ---- Layout ----
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4, 1.4))
  
  op <- par(mar = c(4.6, 4.6, 4.2, 1.2))
  
  # Panel A
  tmpA <- plot_pca(uq_norm_log2, "Before RUV", sex_vec, ageclass_fac)
  mtext("A", side = 3, line = 1.2, adj = 0, font = 2, cex = tmpA$cex_A_B)
  
  # Panel B
  tmpB <- plot_pca(uq_ruv_log2, "After RUV", sex_vec, ageclass_fac)
  mtext("B", side = 3, line = 1.2, adj = 0, font = 2, cex = tmpB$cex_A_B)
  
  # ---- Legends panel ----
  par(mar = c(0, 0, 2.2, 0))
  plot.new()
  
  col_m <- tmpA$col_m
  col_f <- tmpA$col_f
  
  legend(
    x = 0.02, y = 0.6, xpd = NA, bty = "n",
    legend = c("male", "female"),
    pt.bg  = c(col_m, col_f),
    pch    = 21,
    col    = "black",
    pt.cex = 1.8,   
    cex    = 1.25,  
    horiz  = TRUE
  )
  
  K <- length(levels(ageclass_fac))
  legend(
    x = 0.34, y = 0.6, xpd = NA, bty = "n",
    legend = levels(ageclass_fac),
    pch    = 21:(20 + K),
    col    = "black",
    pt.bg  = "grey85",
    pt.cex = 1.8,   
    cex    = 1.15,  
    horiz  = TRUE,
    title  = "Age class",
    title.cex = 1.2 
  )
  
  par(op)
  layout(1)
  
}

save_base_plot(
  draw_pca_ruv,
  file.path(fig_main_dir, "plot_pca_ruv"),
  width_px = 1315, height_px = 789
)

draw_flowchart <- function() {
  op <- par(mar = c(0.5, 0.5, 1.2, 0.5))
  on.exit(par(op))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  boxes <- c(
    "117 participants\n488 quantified proteins",
    "Sample-level QC\nand outlier assessment",
    "Missingness assessment\nand protein filtering",
    "Missing-value imputation\nwithin sex",
    "Upper-quartile normalization\n+ log2 transformation",
    "RUV correction\nwithin sex",
    "117 participants\n370 proteins"
  )
  ys <- seq(0.90, 0.10, length.out = length(boxes))
  for (i in seq_along(boxes)) {
    rect(0.28, ys[i]-0.045, 0.72, ys[i]+0.045, border = "black", lwd = 1.2)
    text(0.50, ys[i], boxes[i], cex = 0.95)
    if (i < length(boxes)) {
      arrows(0.50, ys[i]-0.047, 0.50, ys[i+1]+0.052,
             length = 0.08, lwd = 1.1)
    }
  }
  title("Plasma proteomics preprocessing and quality-control workflow",
        cex.main = 1.2, font.main = 2)
}
save_base_plot(
  draw_flowchart,
  file.path(fig_supp_dir, "flowchart_preprocessing_workflow"),
  width_px = 1400, height_px = 900
)

# Final clean expression matrix and aligned metadata
expr_final <- uq_ruv_log2
meta_final <- meta_aligned

# Check dimensions and manuscript-realized tuning choices.
stopifnot(
  identical(colnames(expr_final), meta_final$Sample_ID),
  nrow(expr_final) == 370L,
  ncol(expr_final) == 117L,
  k_ruv == 8L,
  length(ctl_idx) == 100L
)


reference_rds <- Sys.getenv("NAKO_PREPROCESSING_REFERENCE_RDS", unset = "")
if (nzchar(reference_rds)) {
  if (!file.exists(reference_rds)) {
    stop("NAKO_PREPROCESSING_REFERENCE_RDS does not exist: ", reference_rds)
  }
  ref_expr <- readRDS(reference_rds)
  if (!identical(dim(expr_final), dim(ref_expr)) ||
      !identical(rownames(expr_final), rownames(ref_expr)) ||
      !identical(colnames(expr_final), colnames(ref_expr))) {
    stop("Historical reference and reconstructed matrix have different dimensions/dimnames.")
  }
  delta <- as.numeric(expr_final) - as.numeric(ref_expr)
  cat(sprintf(
    "Reference comparison: max|delta|=%.17g; mean|delta|=%.17g; identical=%s\n",
    max(abs(delta)), mean(abs(delta)), identical(expr_final, ref_expr)
  ))
  if (!isTRUE(all.equal(expr_final, ref_expr, tolerance = 0, check.attributes = TRUE))) {
    stop(
      "The reconstructed matrix is not numerically identical to the supplied historical reference. ",
      "Inspect sessionInfo_preprocessing.txt and the missForest/PhosR versions before downstream analysis."
    )
  }
}

# Save both as RDS (compact and preserves structure)
saveRDS(expr_final, here("analysis", "data", "proteomic_cleaned_expr.rds"))
saveRDS(meta_final, here("analysis", "data", "proteomic_nako_metadata_aligned.rds"))

cat("Saved cleaned expression matrix and metadata:\n",
    " - proteomic_cleaned_expr.rds\n",
    " - proteomic_nako_metadata_aligned.rds\n")

writeLines(capture.output(sessionInfo()),
           here("analysis", "results", "sessionInfo_preprocessing.txt"))

writeLines(
  c(
    sprintf("Final samples: %d", ncol(expr_final)),
    sprintf("Final proteins: %d", nrow(expr_final)),
    sprintf("RUV k: %d", k_ruv),
    sprintf("RUV controls: %d", length(ctl_idx)),
    sprintf("PZP silhouette: %.3f", pzp_sil),
    sprintf("Post-RUV PC1: %.1f%%", ve_post[1]),
    sprintf("Post-RUV PC2: %.1f%%", ve_post[2])
  ),
  here("analysis", "results", "preprocessing_summary.txt")
)

if (nzchar(reference_rds)) {
  message("Preprocessing complete: exact historical-reference check passed.")
} else {
  message("Preprocessing complete: manuscript dimensions/tuning checks passed. Set NAKO_PREPROCESSING_REFERENCE_RDS for exact matrix validation.")
}
