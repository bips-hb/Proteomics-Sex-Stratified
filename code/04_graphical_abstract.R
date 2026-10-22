nako_graphical_abstract <- local({
  PROJECT_ROOT <- NULL
  OUTPUT_BASENAME <- "graphical_abstract_v27"
  FONT_FAMILY <- "Arial"
  REQUIRE_EXACT_FONT <- TRUE
  DPI <- 600
  INITIAL_PROTEINS <- 488L
  INITIAL_PER_STRATUM <- 12L
  EXPORT_SVG <- TRUE
  EXPORT_TOC <- TRUE
  # TOC: the complete main figure, uniformly scaled onto a 50 x 60 mm page.
  # This dependency is used only for the small version; main plotting is unchanged.
  if (EXPORT_TOC && !requireNamespace("grImport2", quietly = TRUE)) {
    stop("For the full-figure TOC, first run install.packages('grImport2').")
  }
  SHOW_IN_PLOT_PANE <- FALSE
  TITLE <- "Sex-specific plasma proteomic signatures associated with chronological aging in the NAKO-Bremen cohort"
  TAKE_HOME <- "Sex-stratified analyses identify distinct age-associated proteomic signatures and biological patterns"
  COL <- list(
    ink = "#123C54", muted = "#45687D", rule = "#D6E0E5", pale = "#F5F8FA",
    male = "#009B6B", female = "#EA9B00", female_text = "#A96C00",
    male_fill = "#DBF1E9", female_fill = "#FCEBCF", shared = "#2D618E", shared_fill = "#E2EEF7",
    Gray_fill = "#EEF2F5", Gray_text = "#7D8B95",
    hallmark_active = "#B7D8F0", hallmark_icon = "#28658D"
  )
  
  tint <- function(color, strength = 0.10) {
    if (length(strength) != 1L || !is.finite(strength) || strength < 0 || strength > 1) {
      stop("Tint strength must be between zero and one.")
    }
    rgb <- grDevices::col2rgb(color)[, 1L] / 255
    rgb <- (1 - strength) + strength * rgb
    grDevices::rgb(rgb[1L], rgb[2L], rgb[3L])
  }
  COL$male_tint <- tint(COL$male, 0.14)
  COL$female_tint <- tint(COL$female, 0.14)
  COL$shared_tint <- tint(COL$shared, 0.12)
  
  script_files <- unlist(lapply(sys.frames(), function(f) if (!is.null(f$ofile)) as.character(f$ofile) else NULL), use.names = FALSE)
  cli_files <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  starts <- unique(c(if (length(script_files)) dirname(script_files), if (length(cli_files)) dirname(cli_files), getwd()))
  if (requireNamespace("here", quietly = TRUE)) starts <- c(starts, here::here())
  is_project <- function(p) all(file.exists(file.path(p, "analysis", "data", c("proteomic_cleaned_expr.rds", "proteomic_nako_metadata_aligned.rds"))))
  find_root <- function() {
    if (!is.null(PROJECT_ROOT)) {
      if (!is_project(PROJECT_ROOT)) stop("PROJECT_ROOT does not contain the required outputs.")
      return(normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE))
    }
    for (s in starts) {
      if (!dir.exists(s)) next
      p <- normalizePath(s, winslash = "/", mustWork = TRUE)
      for (i in seq_len(12L)) {
        if (is_project(p)) return(p)
        q <- dirname(p)
        if (identical(q, p)) break
        p <- q
      }
    }
    stop("Project not found. Open the NAKO project in RStudio or set PROJECT_ROOT.")
  }
  
  root <- find_root()
  data_dir <- file.path(root, "analysis", "data")
  results_dir <- file.path(root, "analysis", "results")
  sets_dir <- file.path(results_dir, "protein_sets")
  output_dir <- file.path(results_dir, "graphical_abstract")
  inputs <- c(
    expression = file.path(data_dir, "proteomic_cleaned_expr.rds"),
    metadata = file.path(data_dir, "proteomic_nako_metadata_aligned.rds"),
    male_set = file.path(sets_dir, "stable_proteins_male_ids.rds"),
    female_set = file.path(sets_dir, "stable_proteins_female_ids.rds")
  )
  if (any(!file.exists(inputs))) stop("Missing required outputs: ", paste(basename(inputs[!file.exists(inputs)]), collapse = ", "))
  
  expr <- as.matrix(readRDS(inputs[["expression"]]))
  meta <- as.data.frame(readRDS(inputs[["metadata"]]))
  if (!is.numeric(expr) || is.null(rownames(expr)) || is.null(colnames(expr))) stop("Invalid expression matrix.")
  if (anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) stop("Duplicate expression identifiers.")
  if (any(!is.finite(expr))) stop("The expression matrix contains non-finite values.")
  id_col <- intersect(c("Sample_ID", "sample_id"), names(meta))[1L]
  age_col <- intersect(c("age_years", "age", "Age", "AGE"), names(meta))[1L]
  if (is.na(id_col) || is.na(age_col) || !"sex" %in% names(meta)) stop("Metadata must contain sample ID, sex and age.")
  if (anyNA(meta[[id_col]]) || anyDuplicated(as.character(meta[[id_col]]))) stop("Missing or duplicate metadata sample identifiers.")
  idx <- match(colnames(expr), as.character(meta[[id_col]]))
  if (anyNA(idx)) stop("Expression samples do not match metadata.")
  meta <- meta[idx, , drop = FALSE]
  sex_raw <- tolower(trimws(as.character(meta$sex)))
  sex <- ifelse(sex_raw %in% c("1", "male", "m", "man", "men"), "male", ifelse(sex_raw %in% c("2", "female", "f", "woman", "women"), "female", NA_character_))
  age <- suppressWarnings(as.numeric(as.character(meta[[age_col]])))
  if (anyNA(sex) || any(!is.finite(age))) stop("Invalid sex or age values.")
  
  get_set <- function(path) {
    x <- as.character(readRDS(path))
    if (anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) stop("Invalid protein set.")
    if (length(setdiff(x, rownames(expr)))) stop("Selected proteins missing from expression matrix.")
    x
  }
  sets <- list(male = get_set(inputs[["male_set"]]), female = get_set(inputs[["female_set"]]))
  shared <- intersect(sets$male, sets$female)
  male_only <- setdiff(sets$male, shared)
  female_only <- setdiff(sets$female, shared)
  
  symbol_map <- c(
    IBP2 = "IGFBP2", LG3BP = "LGALS3BP", PAFA = "PLA2G7", LV469 = "IGLV4-69", TENX = "TNXB",
    ATS13 = "ADAMTS13", IGA2 = "IGHA2", HV307 = "IGHV3-7", CO6A1 = "COL6A1", TRFL = "LTF",
    LV211 = "IGLV2-11", LYSC = "LYZ", IGJ = "JCHAIN", PROP = "CFP", ALS = "IGFALS", GPIX = "GP9",
    ITB3 = "ITGB3", CERU = "CP", ATPA = "ATP5F1A", HV315 = "IGHV3-15", ANGT = "AGT"
  )
  symbols <- function(x) {
    y <- sub("_HUMAN$", "", as.character(x))
    z <- y %in% names(symbol_map)
    y[z] <- unname(symbol_map[y[z]])
    y
  }
  selected_symbols <- list(male = unique(symbols(sets$male)), female = unique(symbols(sets$female)))
  
  quota_labels <- c("20–29", "30–39", "40–49", "50–59", "60–69")
  age_group <- rep(NA_character_, length(age))
  if ("ageclass" %in% names(meta)) {
    ageclass_raw <- trimws(as.character(meta$ageclass))
    ageclass_num <- suppressWarnings(as.integer(ageclass_raw))
    if (all(is.na(ageclass_raw) | !nzchar(ageclass_raw) | (!is.na(ageclass_num) & ageclass_num %in% 1:5))) {
      ok <- !is.na(ageclass_num) & ageclass_num %in% 1:5
      age_group[ok] <- quota_labels[ageclass_num[ok]]
    } else {
      z <- tolower(gsub("[[:space:]]+", "", ageclass_raw))
      age_group[grepl("20.*29|20.*<30|^20$", z)] <- "20–29"
      age_group[grepl("30.*39|30.*<40|^30$", z)] <- "30–39"
      age_group[grepl("40.*49|40.*<50|^40$", z)] <- "40–49"
      age_group[grepl("50.*59|50.*<60|^50$", z)] <- "50–59"
      age_group[grepl("60|\\>=60|≥60|60\\+|60.*69", z)] <- "60–69"
    }
  }
  fallback <- is.na(age_group)
  age_group[fallback & age >= 20 & age < 30] <- "20–29"
  age_group[fallback & age >= 30 & age < 40] <- "30–39"
  age_group[fallback & age >= 40 & age < 50] <- "40–49"
  age_group[fallback & age >= 50 & age < 60] <- "50–59"
  age_group[fallback & age >= 60] <- "60–69"
  if (anyNA(age_group)) stop("Some participants could not be assigned to the five quota-sampling age groups.")
  age_group <- factor(age_group, levels = quota_labels)
  quota <- table(factor(age_group, levels = quota_labels), factor(sex, levels = c("male", "female")))
  
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required.")
  run_limma <- function(g) {
    take <- which(sex == g)
    X <- expr[, take, drop = FALSE]
    age10 <- age[take] / 10
    fit <- limma::eBayes(limma::lmFit(X, stats::model.matrix(~ age10)))
    p <- as.numeric(fit$p.value[, "age10"])
    if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid limma p-values in the ", g, " group.")
    p <- pmax(p, .Machine$double.xmin)
    data.frame(
      protein = rownames(X),
      protein_symbol = symbols(rownames(X)),
      beta = as.numeric(fit$coefficients[, "age10"]),
      pval = p,
      selected = rownames(X) %in% sets[[g]],
      stringsAsFactors = FALSE
    )
  }
  limma_res <- list(male = run_limma("male"), female = run_limma("female"))
  
  compute_ageprot_r2 <- function(g) {
    take <- which(sex == g)
    X <- expr[sets[[g]], take, drop = FALSE]
    if (!is.matrix(X)) X <- matrix(X, nrow = 1, dimnames = list(sets[[g]], colnames(expr)[take]))
    if (any(apply(X, 1L, stats::sd) <= 0)) stop("A selected protein has zero variance in the ", g, " group.")
    Z <- t(scale(t(X), center = TRUE, scale = TRUE))
    if (nrow(Z) == 1L) Z <- matrix(Z, nrow = 1, dimnames = dimnames(X))
    age10 <- age[take] / 10
    design <- cbind("(Intercept)" = 1, age10 = age10)
    beta <- stats::lm.fit(x = design, y = t(Z))$coefficients["age10", ]
    score <- as.numeric(t(Z) %*% beta)
    unname(summary(stats::lm(score ~ age[take]))$r.squared)
  }
  ageprot_r2 <- c(male = compute_ageprot_r2("male"), female = compute_ageprot_r2("female"))
  
  pick_top_proteins <- function(g, hallmark = NULL, n = 4L) {
    d <- limma_res[[g]]
    d <- d[d$selected, , drop = FALSE]
    if (!is.null(hallmark)) {
      keep <- hallmark_map$protein_symbol[hallmark_map$sex == g & hallmark_map$hallmark == hallmark]
      d <- d[d$protein_symbol %in% keep, , drop = FALSE]
    }
    if (!nrow(d) || n <= 0L) return(character(0))
    d <- d[order(d$pval, -abs(d$beta), d$protein_symbol), , drop = FALSE]
    unique(d$protein_symbol)[seq_len(min(n, length(unique(d$protein_symbol))))]
  }
  top_volcano_labels <- list(
    male = pick_top_proteins("male", n = 4L),
    female = pick_top_proteins("female", n = 4L)
  )
  
  hallmark_map <- data.frame(
    sex = c(rep("male", 17), rep("female", 15)),
    protein_symbol = c(
      "IGFBP2", "LGALS3BP", "PIGR", "VASP", "PLA2G7", "IGLV4-69", "LDHA", "SAA4", "TNXB", "PI16", "ADAMTS13", "CD5L", "IGHA2", "IGHV3-7", "COL6A1", "LTF", "IGLV2-11",
      "IGFBP2", "INHBC", "LYZ", "JCHAIN", "CFP", "LDHA", "IGFALS", "GP9", "PRDX6", "ITGB3", "CP", "ATP5F1A", "IGHV3-15", "AGT", "PZP"
    ),
    hallmark = c(
      "Deregulated nutrient sensing", "Chronic inflammation", "Chronic inflammation", "Loss of proteostasis", "Altered intercellular communication", "Chronic inflammation", "Deregulated nutrient sensing", "Chronic inflammation", "Loss of proteostasis", "Loss of proteostasis", "Loss of proteostasis", "Chronic inflammation", "Chronic inflammation", "Chronic inflammation", "Loss of proteostasis", "Mitochondrial dysfunction", "Chronic inflammation",
      "Deregulated nutrient sensing", "Deregulated nutrient sensing", "Chronic inflammation", "Loss of proteostasis", "Deregulated nutrient sensing", "Mitochondrial dysfunction", "Deregulated nutrient sensing", "Altered intercellular communication", "Mitochondrial dysfunction", "Loss of proteostasis", "Mitochondrial dysfunction", "Mitochondrial dysfunction", "Chronic inflammation", "Altered intercellular communication", "Loss of proteostasis"
    ),
    stringsAsFactors = FALSE
  )
  
  hallmarks12 <- data.frame(
    hallmark = c(
      "Genomic instability", "Telomere attrition", "Epigenetic alterations", "Loss of proteostasis", "Disabled macroautophagy",
      "Deregulated nutrient sensing", "Mitochondrial dysfunction", "Cellular senescence", "Stem cell exhaustion",
      "Altered intercellular communication", "Chronic inflammation", "Dysbiosis"
    ),
    label = c(
      "Genomic\ninstability", "Telomere\nattrition", "Epigenetic\nalterations", "Loss of\nproteostasis", "Disabled\nmacroautophagy",
      "Deregulated\nnutrient sensing", "Mitochondrial\ndysfunction", "Cellular\nsenescence", "Stem cell\nexhaustion",
      "Altered\nintercellular\ncommunication", "Chronic\ninflammation", "Dysbiosis"
    ),
    stringsAsFactors = FALSE
  )
  
  represented_hallmarks <- c("Deregulated nutrient sensing", "Chronic inflammation", "Loss of proteostasis", "Altered intercellular communication", "Mitochondrial dysfunction")
  hallmarks12 <- hallmarks12[
    order(!(hallmarks12$hallmark %in% represented_hallmarks)), , drop = FALSE
  ]
  hallmarks12$fill <- ifelse(
    hallmarks12$hallmark %in% represented_hallmarks, COL$hallmark_active, COL$Gray_fill
  )
  hallmark_protein_lines <- setNames(vector("list", nrow(hallmarks12)), hallmarks12$hallmark)
  for (h in hallmarks12$hallmark) {
    male_p <- pick_top_proteins("male", h, if (h %in% represented_hallmarks) 2L else 0L)
    female_p <- pick_top_proteins("female", h, if (h %in% represented_hallmarks) 2L else 0L)
    if (h == "Deregulated nutrient sensing") {
      if (!"IGFBP2" %in% male_p && "IGFBP2" %in% selected_symbols$male) male_p <- c("IGFBP2", male_p)
      if ("INHBC" %in% selected_symbols$female) female_p <- unique(c("INHBC", female_p))
    }
    hallmark_protein_lines[[h]] <- list(
      male = head(unique(male_p[nzchar(male_p)]), 2L),
      female = head(unique(female_p[nzchar(female_p)]), 2L)
    )
  }
  
  n_male <- sum(sex == "male")
  n_female <- sum(sex == "female")
  initial_quota <- matrix(
    INITIAL_PER_STRATUM, nrow = length(quota_labels), ncol = 2L,
    dimnames = list(quota_labels, c("male", "female"))
  )
  initial_participants <- sum(initial_quota)
  initial_n_male <- sum(initial_quota[, "male"])
  initial_n_female <- sum(initial_quota[, "female"])
  excluded_after_qc <- initial_participants - ncol(expr)
  if (any(quota > initial_quota) || excluded_after_qc < 0) {
    stop("The analytical sample exceeds the initial quota design.")
  }
  venn_box_proteins <- lapply(selected_symbols, sort)
  for (g in names(venn_box_proteins)) {
    if (length(venn_box_proteins[[g]]) != length(sets[[g]]) ||
        !setequal(venn_box_proteins[[g]], symbols(sets[[g]]))) {
      stop("Venn protein boxes must include every selected protein exactly once.")
    }
  }
  if (!capabilities("cairo")) stop("This script requires R with Cairo graphics support.")
  if (!is.numeric(DPI) || length(DPI) != 1L || !is.finite(DPI) || DPI < 300) stop("DPI must be at least 300.")
  if (!requireNamespace("systemfonts", quietly = TRUE)) stop("Install package 'systemfonts' to verify the requested font.")
  available_fonts <- unique(systemfonts::system_fonts()$family)
  requested_font <- match(tolower(FONT_FAMILY), tolower(available_fonts))
  if (is.na(requested_font)) {
    if (REQUIRE_EXACT_FONT) stop(
      "Font '", FONT_FAMILY, "' was not found. Install it on this computer and restart R. ",
      "No substitute font has been used."
    )
    alternatives <- c("Arial", "Liberation Sans", "Arimo", "Helvetica", "Nimbus Sans", "sans")
    found <- alternatives[tolower(alternatives) %in% tolower(available_fonts)]
    resolved_font <- if (length(found)) found[1L] else "serif"
    warning("Preview font substituted: ", resolved_font, ". This is not Times New Roman.", call. = FALSE)
  } else {
    resolved_font <- available_fonts[requested_font]
  }
  
  measure_text <- function(label, size, face = "plain", lineheight = 1.05) {
    g <- grid::textGrob(label, gp = grid::gpar(
      fontfamily = resolved_font, fontsize = size, fontface = face, lineheight = lineheight
    ))
    c(width = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
      height = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE))
  }
  
  fitted_font_size <- function(label, size, width_mm = Inf, height_mm = Inf,
                               face = "plain", lineheight = 1.05) {
    if (length(size) != 1L || !is.finite(size) || size <= 0) stop("Text size must be a positive number.")
    if (length(width_mm) != 1L || is.na(width_mm) || width_mm <= 0 ||
        length(height_mm) != 1L || is.na(height_mm) || height_mm <= 0) stop("Text dimensions must be positive.")
    current <- size
    for (attempt in seq_len(12L)) {
      measured <- measure_text(label, current, face, lineheight)
      ratios <- c(width_mm / max(measured["width"], 1e-8),
                  height_mm / max(measured["height"], 1e-8))
      factor <- min(1, ratios)
      if (factor >= 0.999) return(current)
      current <- current * factor * 0.985
      if (!is.finite(current) || current <= 0) stop("Invalid text dimensions for: ", label)
    }
    current
  }
  
  draw <- function() {
    W <- 360
    H <- 150
    panel_bottom <- 10
    panel_limits <- rbind(
      study = c(12, 70), pipeline = c(80, 132),
      signatures = c(142, 246), hallmarks = c(256, 348)
    )
    panel_centres <- rowMeans(panel_limits)
    panel_widths <- panel_limits[, 2L] - panel_limits[, 1L]
    separators <- (panel_limits[1:3, 2L] + panel_limits[2:4, 1L]) / 2
    Y <- function(y) y * H / 150
    u <- function(x) grid::unit(x, "native")
    text_at <- function(label, x, y, size = 9, color = COL$ink, face = "plain",
                        just = "centre", rot = 0, width = NULL, height = NULL, lineheight = 1.05) {
      width_mm <- if (is.null(width)) Inf else grid::convertWidth(u(width), "mm", valueOnly = TRUE)
      height_mm <- if (is.null(height)) Inf else grid::convertHeight(u(Y(height)), "mm", valueOnly = TRUE)
      size <- fitted_font_size(label, size, width_mm, height_mm, face, lineheight)
      grid::grid.text(label, u(x), u(Y(y)), just = just, rot = rot, gp = grid::gpar(
        fontfamily = resolved_font, fontsize = size, col = color, fontface = face, lineheight = lineheight
      ))
    }
    
    line <- function(x, y, color = COL$rule, lwd = 1, arrow = FALSE) {
      grid::grid.lines(u(x), u(Y(y)), gp = grid::gpar(col = color, lwd = lwd, lineend = "round"), arrow = if (arrow) grid::arrow(length = grid::unit(1.5, "mm"), type = "closed") else NULL)
    }
    panel_arrow <- function(mid, y = 76) {
      arrow_length <- 8.4
      arrow_height <- 6.8
      shaft_height <- 2.6
      head_length <- 3.3
      x1 <- mid - arrow_length / 2
      x2 <- mid + arrow_length / 2
      shoulder <- x2 - head_length
      grid::grid.polygon(
        u(c(x1, shoulder, shoulder, x2, shoulder, shoulder, x1)),
        u(Y(y + c(shaft_height / 2, shaft_height / 2, arrow_height / 2, 0,
                  -arrow_height / 2, -shaft_height / 2, -shaft_height / 2))),
        gp = grid::gpar(fill = COL$muted, col = NA)
      )
    }
    box <- function(x, y, w, h, fill = "white", stroke = COL$rule, lwd = 1, radius = 1) {
      grid::grid.roundrect(u(x), u(Y(y)), u(w), u(Y(h)), r = grid::unit(radius, "mm"), gp = grid::gpar(fill = fill, col = stroke, lwd = lwd))
    }
    circle <- function(x, y, r, fill, stroke = NA, lwd = 1) {
      grid::grid.circle(u(x), u(Y(y)), r = u(r), gp = grid::gpar(fill = fill, col = stroke, lwd = lwd))
    }
    badge <- function(x, y, w, h, fill) {
      grid::grid.rect(u(x), u(Y(y)), u(w), u(Y(h)),
                      gp = grid::gpar(fill = fill, col = NA, lwd = 0))
    }
    medical_silhouette <- function(x, y, group, height = 5.1) {
      if (!group %in% c("male", "female")) stop("Unknown silhouette group.")
      color <- COL[[group]]
      circle(x, y + height * 0.37, height * 0.095, color)
      if (group == "male") {
        right <- rbind(
          c(0.052, 0.245), c(0.174, 0.216), c(0.211, 0.173),
          c(0.257, -0.125), c(0.246, -0.172), c(0.208, -0.181),
          c(0.190, -0.146), c(0.154, 0.056), c(0.143, -0.151),
          c(0.129, -0.485), c(0.043, -0.485), c(0.029, -0.195)
        )
      } else {
        right <- rbind(
          c(0.050, 0.245), c(0.146, 0.210), c(0.181, 0.166),
          c(0.238, -0.122), c(0.229, -0.167), c(0.195, -0.177),
          c(0.177, -0.145), c(0.136, 0.061), c(0.098, -0.054),
          c(0.145, -0.162), c(0.123, -0.485), c(0.041, -0.485),
          c(0.026, -0.206)
        )
      }
      left <- right[nrow(right):1L, , drop = FALSE]
      left[, 1L] <- -left[, 1L]
      outline <- rbind(right, left)
      grid::grid.polygon(u(x + height * outline[, 1L]),
                         u(Y(y + height * outline[, 2L])),
                         gp = grid::gpar(fill = color, col = NA, lwd = 0))
      invisible(NULL)
    }
    venn <- function(cx, cy, r) {
      left <- cx - r * 0.58
      right <- cx + r * 0.58
      circle(left, cy, r, grDevices::adjustcolor(COL$male, alpha.f = 0.22), COL$male, 1)
      circle(right, cy, r, grDevices::adjustcolor(COL$female, alpha.f = 0.22), COL$female, 1)
      text_at(as.character(length(male_only)), left - r * 0.32, cy, 12.2, COL$male, "bold")
      text_at(as.character(length(shared)), cx, cy, 12.2, COL$shared, "bold")
      text_at(as.character(length(female_only)), right + r * 0.32, cy, 12.2, COL$female_text, "bold")
    }
    
    volcano_xr <- lapply(limma_res, function(z) {
      xr <- range(z$beta, finite = TRUE)
      xr + c(-1, 1) * max(diff(xr) * 0.06, 0.01)
    })
    volcano_x_ticks <- lapply(names(limma_res), function(g) {
      ticks <- pretty(volcano_xr[[g]], n = 5)
      ticks[ticks >= volcano_xr[[g]][1] & ticks <= volcano_xr[[g]][2]]
    })
    names(volcano_x_ticks) <- names(limma_res)
    volcano_yr <- lapply(limma_res, function(z) {
      top <- max(-log10(z$pval), finite = TRUE) * 1.05
      c(0, top)
    })
    volcano_y_ticks <- lapply(names(limma_res), function(g) {
      ticks <- pretty(volcano_yr[[g]], n = 4)
      ticks[ticks >= volcano_yr[[g]][1] & ticks <= volcano_yr[[g]][2]]
    })
    names(volcano_y_ticks) <- names(limma_res)
    
    volcano <- function(g, cx, cy, w, h, show_y_ticks = TRUE) {
      d <- limma_res[[g]]
      sel <- d[d$selected, , drop = FALSE]
      color <- if (g == "male") COL$male else COL$female
      label_color <- if (g == "male") COL$male else COL$female_text
      label_targets <- top_volcano_labels[[g]]
      x_range <- volcano_xr[[g]]
      x_ticks <- volcano_x_ticks[[g]]
      y_range <- volcano_yr[[g]]
      y_ticks <- volcano_y_ticks[[g]]
      
      grid::pushViewport(grid::viewport(x = u(cx), y = u(Y(cy)), width = u(w), height = u(Y(h)), xscale = x_range, yscale = y_range, clip = "on"))
      grid::grid.points(u(d$beta), u(-log10(d$pval)), pch = 16, size = grid::unit(0.68, "mm"), gp = grid::gpar(col = grDevices::adjustcolor(COL$muted, 0.27)))
      grid::grid.points(u(sel$beta), u(-log10(sel$pval)), pch = 21, size = grid::unit(1.15, "mm"), gp = grid::gpar(col = color, fill = grDevices::adjustcolor(color, 0.84), lwd = 0.7))
      grid::popViewport()
      
      line(c(cx - w / 2, cx - w / 2), c(cy - h / 2, cy + h / 2), COL$muted, 0.7)
      line(c(cx - w / 2, cx + w / 2), c(cy - h / 2, cy - h / 2), COL$muted, 0.7)
      line(c(cx + w / 2, cx + w / 2), c(cy - h / 2, cy + h / 2), COL$rule, 0.5)
      line(c(cx - w / 2, cx + w / 2), c(cy + h / 2, cy + h / 2), COL$rule, 0.5)
      
      for (v in x_ticks) {
        xx <- cx - w / 2 + (v - x_range[1]) / diff(x_range) * w
        line(c(xx, xx), c(cy - h / 2, cy - h / 2 - 1.2), COL$muted, 0.6)
        text_at(format(round(v, 1), nsmall = 1, trim = TRUE), xx, cy - h / 2 - 3.8, 7.2, COL$muted)
      }
      if (show_y_ticks) {
        for (v in y_ticks) {
          yy <- cy - h / 2 + (v - y_range[1]) / diff(y_range) * h
          line(c(cx - w / 2, cx - w / 2 - 1.2), c(yy, yy), COL$muted, 0.6)
          text_at(format(round(v, 1), nsmall = 1, trim = TRUE), cx - w / 2 - 3.1, yy, 7.2, COL$muted, just = "right")
        }
      }
      
      lab_df <- sel[match(label_targets, sel$protein_symbol), , drop = FALSE]
      lab_df <- lab_df[!is.na(lab_df$protein_symbol), , drop = FALSE]
      placed <- list()
      point_x <- cx - w / 2 + (d$beta - x_range[1]) / diff(x_range) * w
      point_y <- cy - h / 2 + (-log10(d$pval) - y_range[1]) / diff(y_range) * h
      disjoint <- function(a, b, gap = 0.6) {
        a[2] + gap <= b[1] || b[2] + gap <= a[1] || a[4] + gap <= b[3] || b[4] + gap <= a[3]
      }
      if (nrow(lab_df)) for (i in seq_len(nrow(lab_df))) {
        xx <- cx - w / 2 + (lab_df$beta[i] - x_range[1]) / diff(x_range) * w
        yy <- cy - h / 2 + (-log10(lab_df$pval[i]) - y_range[1]) / diff(y_range) * h
        label <- lab_df$protein_symbol[i]
        dims <- measure_text(label, 7.0, "bold")
        lw <- unname(dims["width"]) + 1.10
        lh <- unname(dims["height"]) + 0.55
        if (g == "male" && label == "VASP") {
          candidates <- data.frame(
            dx = c(-(lw / 2 + 2.0), -(lw / 2 + 4.0), -(lw / 2 + 6.0), -(lw / 2 + 8.0),
                   -(lw / 2 + 2.0), -(lw / 2 + 4.0), -(lw / 2 + 6.0),
                   (lw / 2 + 2.0),  (lw / 2 + 4.0)),
            dy = c(0.0, 0.0, 0.0, 0.0, 2.0, -2.0, 3.6, 0.0, 2.0)
          )
        } else {
          candidates <- expand.grid(dx = c(-(lw / 2 + c(1.0, 2.6, 5.0)), lw / 2 + c(1.0, 2.6, 5.0)),
                                    dy = c(2.0, -2.0, 3.6, -3.6, 5.2, -5.2, 7.0, -7.0, 0))
        }
        candidates$x <- pmax(cx - w / 2 + lw / 2 + 0.5,
                             pmin(cx + w / 2 - lw / 2 - 0.5, xx + candidates$dx))
        candidates$y <- pmax(cy - h / 2 + lh / 2 + 0.5,
                             pmin(cy + h / 2 - lh / 2 - 0.5, yy + candidates$dy))
        bounds <- lapply(seq_len(nrow(candidates)), function(j) c(
          candidates$x[j] - lw / 2, candidates$x[j] + lw / 2,
          candidates$y[j] - lh / 2, candidates$y[j] + lh / 2
        ))
        costs <- vapply(seq_along(bounds), function(j) {
          b <- bounds[[j]]
          if (length(placed) && !all(vapply(placed, function(old) disjoint(b, old), logical(1)))) return(Inf)
          covers <- point_x > b[1] & point_x < b[2] & point_y > b[3] & point_y < b[4]
          pref_penalty <- if (g == "male" && label == "VASP" && candidates$x[j] > xx) 1000 else 0
          (candidates$x[j] - xx)^2 + (candidates$y[j] - yy)^2 + pref_penalty +
            20 * sum(covers) + 120 * sum(covers & d$selected)
        }, numeric(1))
        if (!any(is.finite(costs))) next
        j <- which.min(costs)
        bx <- candidates$x[j]
        by <- candidates$y[j]
        b <- bounds[[j]]
        end_x <- max(b[1], min(xx, b[2]))
        end_y <- max(b[3], min(yy, b[4]))
        line(c(xx, end_x), c(yy, end_y), label_color, 0.35)
        text_at(label, bx, by, 7.0, label_color, "bold")
        placed[[length(placed) + 1L]] <- b
      }
    }
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(xscale = c(0, W), yscale = c(0, H), gp = grid::gpar(fontfamily = resolved_font)))
    on.exit(grid::popViewport(), add = TRUE)
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    
    for (i in seq_len(nrow(panel_limits))) {
      heading <- c("1  Study design", "2  Analytical pipeline",
                   "3  Age-associated signatures", "4  Biological interpretation")[i]
      text_at(heading, panel_centres[i], 142, 12.0, COL$ink, "bold",
              width = panel_widths[i] - 2)
    }
    for (mid in separators) {
      line(rep(mid, 2), c(4, 146), COL$rule, 0.8)
      panel_arrow(mid, y = 80)
    }
    
    study_cx <- panel_centres[["study"]]
    age_x <- 22
    quota_x <- c(male = 44, female = 62)
    text_at("NAKO-Bremen", study_cx, 131, 12.5, COL$ink, "bold")
    text_at(sprintf("%d participants", initial_participants), study_cx, 124, 10.8, COL$muted)
    text_at("Age group\n(years)", age_x, 112, 9.6, COL$muted, "bold", width = 20)
    text_at("Male", quota_x["male"], 112, 9.6, COL$male, "bold", width = 16)
    text_at("Female", quota_x["female"], 112, 9.6, COL$female_text, "bold", width = 16)
    ys <- c(96, 80, 64, 48, 32)
    for (i in seq_along(ys)) {
      text_at(rownames(initial_quota)[i], age_x, ys[i], 9.5, COL$muted)
      for (g in names(quota_x)) {
        fill <- if (g == "male") COL$male_tint else COL$female_tint
        badge(quota_x[g], ys[i], 16, 6.2, fill)
        medical_silhouette(quota_x[g] - 4.5, ys[i], g)
        text_at(paste0("n=", initial_quota[i, g]), quota_x[g] + 2, ys[i],
                9.2, COL$ink, "bold")
      }
    }
    total_y <- 16
    text_at("Total", age_x, total_y, 9.2, COL$muted, "bold", c("centre", "bottom"))
    text_at(sprintf("n=%d", initial_n_male), quota_x["male"], total_y,
            9.5, COL$male, "bold", c("centre", "bottom"))
    text_at(sprintf("n=%d", initial_n_female), quota_x["female"], total_y,
            9.5, COL$female_text, "bold", c("centre", "bottom"))
    
    pipeline_cx <- panel_centres[["pipeline"]]
    pipe_y <- c(122, 101, 80, 59, 38, 17)
    pipe_lab <- c(
      sprintf("LC-MS/MS profiling\n%d proteins quantified\n%d participants retained after QC", INITIAL_PROTEINS, ncol(expr)),
      sprintf("Quality control, imputation\nand normalization\n%d proteins retained", nrow(expr)),
      "Sex-stratified ordinal elastic net\nwith nested cross-validation",
      "Unsupervised internal assessment\nof held-out samples",
      "Stable age-associated\nprotein signatures",
      "AgeProt score\nWeighted sum of selected proteins"
    )
    pipe_font_size <- min(vapply(pipe_lab, function(label) {
      fitted_font_size(label, 9.0, width_mm = panel_widths["pipeline"] - 4)
    }, numeric(1)))
    pipe_sizes <- rep(pipe_font_size, length(pipe_lab))
    pipe_heights <- vapply(seq_along(pipe_lab), function(i) {
      unname(measure_text(pipe_lab[i], pipe_sizes[i])["height"])
    }, numeric(1))
    pipe_y[length(pipe_y)] <- 17
    for (i in seq_along(pipe_y)) {
      text_at(pipe_lab[i], pipeline_cx, pipe_y[i], pipe_sizes[i], COL$muted,
              width = panel_widths["pipeline"] - 4)
      if (i < length(pipe_y)) {
        arrow_top <- pipe_y[i] - pipe_heights[i] / 2 - 1.4
        arrow_bottom <- pipe_y[i + 1L] + pipe_heights[i + 1L] / 2 + 1.4
        if (arrow_top <= arrow_bottom) stop("Pipeline steps overlap.")
        line(rep(pipeline_cx, 2), c(arrow_top, arrow_bottom), COL$muted, 1, TRUE)
      }
    }
    
    signatures_cx <- panel_centres[["signatures"]]
    volcano_width <- 40.4
    volcano_height <- 33
    volcano_y <- 96
    volcano_left <- panel_limits["signatures", 1L] + 11.5
    volcano_gap <- 8
    volcano_centres <- c(volcano_left, volcano_left + volcano_width + volcano_gap) + volcano_width / 2
    y_title_x <- panel_limits["signatures", 1L] + 2
    badge(signatures_cx, 128.5, panel_widths["signatures"] - 2, 9.8, COL$shared_tint)
    text_at("AgeProt and chronological age", signatures_cx, 130.6, 9.5, COL$shared, "bold",
            width = panel_widths["signatures"] - 8)
    text_at(sprintf("R² = %.2f", ageprot_r2["male"]), volcano_centres[1], 126.3,
            8.8, COL$male, "bold", width = 38)
    text_at(sprintf("R² = %.2f", ageprot_r2["female"]), volcano_centres[2], 126.3,
            8.8, COL$female_text, "bold", width = 38)
    volcano("male", volcano_centres[1], volcano_y, volcano_width, volcano_height, show_y_ticks = TRUE)
    volcano("female", volcano_centres[2], volcano_y, volcano_width, volcano_height, show_y_ticks = TRUE)
    text_at("Age association (β per 10 years)", mean(volcano_centres), 69.2, 8.2, COL$muted)
    text_at("−log10(p-value)", y_title_x, volcano_y, 8.0, COL$muted, rot = 90)
    
    venn_cx <- signatures_cx
    venn_cy <- 34
    venn_radius <- 8.5
    protein_box_width <- 33
    protein_box_height <- 33
    protein_box_x <- c(male = signatures_cx - 35, female = signatures_cx + 35)
    protein_column_width <- (protein_box_width - 6) / 2
    protein_columns <- lapply(venn_box_proteins, function(proteins) {
      n_rows <- ceiling(length(proteins) / 2)
      split(proteins, ceiling(seq_along(proteins) / n_rows))
    })
    protein_box_font <- min(vapply(unlist(protein_columns, recursive = FALSE), function(proteins) {
      fitted_font_size(paste(proteins, collapse = "\n"), 7.6, width_mm = protein_column_width,
                       height_mm = 25, face = "bold", lineheight = 1.06)
    }, numeric(1)))
    for (g in names(protein_box_x)) {
      color <- if (g == "male") COL$male else COL$female_text
      fill <- if (g == "male") COL$male_tint else COL$female_tint
      text_at(if (g == "male") "Male" else "Female", protein_box_x[g],
              venn_cy + protein_box_height / 2 + 4.5, 9.6, color, "bold",
              width = protein_box_width - 4)
      box(protein_box_x[g], venn_cy, protein_box_width, protein_box_height,
          fill = fill, stroke = NA, lwd = 0, radius = 1.2)
      text_at(sprintf("%d proteins", length(venn_box_proteins[[g]])), protein_box_x[g],
              venn_cy + 13.2, 8.5, color, "bold", width = protein_box_width - 4)
      for (j in seq_along(protein_columns[[g]])) {
        column_x <- protein_box_x[g] + (j - 1.5) * (protein_column_width + 2)
        text_at(paste(protein_columns[[g]][[j]], collapse = "\n"), column_x, venn_cy + 10,
                protein_box_font, color, "bold", just = c("centre", "top"),
                width = protein_column_width, height = 25, lineheight = 1.06)
      }
    }
    venn(venn_cx, venn_cy, venn_radius)
    line(c(protein_box_x["male"] + protein_box_width / 2 + 1,
           venn_cx - 1.58 * venn_radius - 0.7), rep(venn_cy, 2), COL$male, 0.8, TRUE)
    line(c(protein_box_x["female"] - protein_box_width / 2 - 1,
           venn_cx + 1.58 * venn_radius + 0.7), rep(venn_cy, 2), COL$female, 0.8, TRUE)
    text_at(paste0("Shared proteins: ", paste(symbols(shared), collapse = " and ")),
            signatures_cx, panel_bottom, 8.9, COL$shared, "bold", just = c("centre", "bottom"), width = 98)
    
    cx <- panel_centres[["hallmarks"]] + 1.5
    cy <- 72
    r1 <- 49
    hallmark_scale <- r1 / 38
    r0 <- 10 * hallmark_scale
    angs <- seq(pi / 2, pi / 2 - 2 * pi, length.out = nrow(hallmarks12) + 1L)
    icon_radius <- 18.8
    icon_size <- 7.0
    text_inner_radius <- icon_radius + icon_size / 2 + 0.7
    draw_hallmark_icon <- function(hallmark, x, y, size, angle, active) {
      icon_ink <- if (active) COL$hallmark_icon else COL$Gray_text
      icon_fill <- tint(icon_ink, 0.14)
      grid::pushViewport(grid::viewport(
        x = u(x), y = u(Y(y)), width = grid::unit(size, "mm"),
        height = grid::unit(size, "mm"), xscale = c(-1, 1), yscale = c(-1, 1),
        angle = angle, clip = "off"
      ))
      on.exit(grid::popViewport(), add = TRUE)
      iu <- function(v) grid::unit(v, "native")
      path <- function(px, py, width = 1.0, color = icon_ink, arrow = FALSE) {
        grid::grid.lines(iu(px), iu(py), gp = grid::gpar(
          col = color, lwd = width, lineend = "round", linejoin = "round"
        ), arrow = if (arrow) grid::arrow(
          length = grid::unit(0.65, "mm"), type = "closed"
        ) else NULL)
      }
      polygon <- function(px, py, fill = icon_fill, width = 1.0) {
        grid::grid.polygon(iu(px), iu(py), gp = grid::gpar(
          fill = fill, col = icon_ink, lwd = width, linejoin = "round"
        ))
      }
      ellipse <- function(ex, ey, rx, ry = rx, fill = icon_fill, width = 1.0) {
        a <- seq(0, 2 * pi, length.out = 80)
        polygon(ex + rx * cos(a), ey + ry * sin(a), fill, width)
      }
      dot <- function(dx, dy, radius = 0.07) {
        ellipse(dx, dy, radius, fill = icon_ink, width = 0.35)
      }
      if (hallmark == "Genomic instability") {
        t <- seq(-0.82, 0.82, length.out = 140)
        z <- 0.40 * sin((t + 0.82) * 2 * pi / 1.64)
        path(t, z, 1.3)
        keep <- abs(t) > 0.12
        for (side in c(-1, 1)) {
          take <- keep & sign(t) == side
          path(t[take], -z[take], 1.3)
        }
        for (v in seq(-0.70, 0.70, by = 0.20)) {
          zy <- 0.40 * sin((v + 0.82) * 2 * pi / 1.64)
          path(c(v, v), c(zy, -zy), 0.65)
        }
      } else if (hallmark == "Telomere attrition") {
        path(c(-0.48, -0.20, 0.20, 0.48), c(0.78, 0.22, -0.22, -0.78), 3.3)
        path(c(0.48, 0.20, -0.20, -0.48), c(0.78, 0.22, -0.22, -0.78), 3.3)
        for (sx in c(-1, 1)) {
          path(sx * c(0.44, 0.49), c(0.70, 0.80), 3.5, icon_fill)
          path(sx * c(0.47, 0.49), c(-0.75, -0.80), 3.5, icon_fill)
        }
      } else if (hallmark == "Epigenetic alterations") {
        for (ex in c(-0.40, 0.40)) {
          ellipse(ex, 0, 0.30, 0.38)
          a <- seq(-pi / 2, 3 * pi / 2, length.out = 100)
          path(ex + 0.33 * cos(a), 0.28 * sin(a) + 0.06, 0.75)
          path(c(ex, ex), c(0.37, 0.66), 0.75)
          dot(ex, 0.71, 0.09)
        }
        path(c(-0.87, -0.70, -0.38, 0.00, 0.38, 0.70, 0.87),
             c(-0.38, -0.16, 0.06, -0.12, 0.08, 0.18, 0.38), 0.95)
      } else if (hallmark == "Loss of proteostasis") {
        ribbon_curve <- function(p0, p1, p2, p3, width = 1.0) {
          t <- seq(0, 1, length.out = 70)
          xy <- outer((1 - t)^3, p0) + outer(3 * (1 - t)^2 * t, p1) +
            outer(3 * (1 - t) * t^2, p2) + outer(t^3, p3)
          path(xy[, 1], xy[, 2], width)
        }
        ribbon_curve(c(-0.84, -0.65), c(-0.83, -0.38),
                     c(-0.71, -0.62), c(-0.54, -0.48))
        t <- seq(0, 1, length.out = 130)
        helix_x <- -0.54 + 0.14 * sin(4 * pi * t)
        helix_y <- -0.48 + 0.96 * t
        path(helix_x, helix_y, 2.5, icon_fill)
        path(helix_x, helix_y, 1.15)
        ribbon_curve(c(-0.54, 0.48), c(-0.20, 0.83),
                     c(-0.30, -0.63), c(0.10, -0.43))
        polygon(c(0.02, 0.18, 0.18, 0.30, 0.10, -0.10, 0.02),
                c(-0.43, -0.43, 0.28, 0.28, 0.55, 0.28, 0.28),
                fill = icon_fill, width = 0.8)
        ribbon_curve(c(0.10, 0.55), c(0.17, 0.83),
                     c(0.57, 0.83), c(0.57, 0.50))
        polygon(c(0.49, 0.65, 0.65, 0.77, 0.57, 0.37, 0.49),
                c(0.50, 0.50, -0.26, -0.26, -0.54, -0.26, -0.26),
                fill = icon_fill, width = 0.8)
        ribbon_curve(c(0.57, -0.54), c(0.59, -0.83),
                     c(0.01, -0.79), c(-0.10, -0.56))
      } else if (hallmark == "Disabled macroautophagy") {
        ellipse(0, 0, 0.78, 0.70, fill = NA)
        ellipse(0, 0, 0.61, 0.54, fill = NA, width = 0.65)
        ellipse(-0.14, -0.05, 0.30, 0.17)
        path(c(-0.34, -0.22, -0.12, -0.02, 0.08),
             c(-0.03, 0.03, -0.12, 0.02, -0.08), 0.65)
        dot(0.30, 0.21, 0.08)
        path(c(0.48, 0.83), c(-0.78, -0.43), 1.5)
        path(c(0.48, 0.83), c(-0.43, -0.78), 1.5)
      } else if (hallmark == "Deregulated nutrient sensing") {
        path(c(-0.44, -0.44), c(-0.78, 0.14), 1.8)
        path(c(-0.68, -0.68, -0.44, -0.20, -0.20),
             c(0.75, 0.20, 0.10, 0.20, 0.75), 1.2)
        path(c(-0.44, -0.44), c(0.75, 0.22), 1.1)
        polygon(c(0.39, 0.50, 0.60, 0.67, 0.69, 0.69, 0.39),
                c(0.80, 0.73, 0.59, 0.38, 0.16, -0.10, -0.10),
                fill = icon_ink, width = 0.55)
        path(c(0.46, 0.46), c(-0.08, -0.78), 2.2)
      } else if (hallmark == "Mitochondrial dysfunction") {
        a <- seq(0, 2 * pi, length.out = 120)
        polygon(0.85 * cos(a), 0.49 * sin(a) * (1 + 0.12 * cos(a)))
        path(c(-0.65, -0.51, -0.42, -0.31, -0.22, -0.11,
               0.01, 0.13, 0.25, 0.37, 0.49, 0.64),
             c(0.08, 0.26, -0.17, -0.27, 0.25, 0.28,
               -0.24, -0.28, 0.25, 0.20, -0.19, -0.03), 1.0)
      } else if (hallmark == "Cellular senescence") {
        a <- seq(0, 2 * pi, length.out = 120)
        radius <- 0.61 + 0.12 * sin(5 * a) + 0.06 * cos(9 * a)
        polygon(radius * cos(a), 0.87 * radius * sin(a))
        ellipse(-0.12, 0.03, 0.21, 0.18, fill = icon_ink, width = 0.35)
        for (q in list(c(0.26, 0.21), c(0.34, -0.13), c(-0.29, -0.27))) {
          dot(q[1], q[2], 0.055)
        }
      } else if (hallmark == "Stem cell exhaustion") {
        cell_centre <- c(-0.22, -0.22)
        cell_radius <- 0.50
        ellipse(cell_centre[1], cell_centre[2], cell_radius)
        ellipse(-0.27, -0.22, 0.22, fill = icon_ink, width = 0.35)
        loop_centre <- c(0.23, 0.28)
        loop_radius <- 0.54
        delta <- cell_centre - loop_centre
        distance <- sqrt(sum(delta^2))
        direction <- delta / distance
        along <- (loop_radius^2 - cell_radius^2 + distance^2) / (2 * distance)
        across <- sqrt(loop_radius^2 - along^2)
        foot <- loop_centre + along * direction
        normal <- c(-direction[2], direction[1])
        p0 <- foot + across * normal
        p3 <- foot - across * normal
        start_angle <- atan2(p0[2] - loop_centre[2], p0[1] - loop_centre[1])
        end_angle <- atan2(p3[2] - loop_centre[2], p3[1] - loop_centre[1])
        if (end_angle < start_angle) end_angle <- end_angle + 2 * pi
        a <- seq(start_angle, end_angle, length.out = 100)
        path(loop_centre[1] + loop_radius * cos(a),
             loop_centre[2] + loop_radius * sin(a), 0.95)
        tangent <- c(-sin(end_angle), cos(end_angle))
        normal <- c(-tangent[2], tangent[1])
        head_base <- p3 - 0.19 * tangent
        head <- rbind(p3, head_base + 0.10 * normal, head_base - 0.10 * normal)
        polygon(head[, 1], head[, 2], fill = icon_ink, width = 0.35)
      } else if (hallmark == "Altered intercellular communication") {
        ellipse(-0.57, -0.09, 0.30, 0.37)
        ellipse(0.57, -0.09, 0.30, 0.37)
        dot(-0.59, -0.08, 0.11)
        dot(0.59, -0.08, 0.11)
        for (dx in c(-0.16, 0, 0.16)) dot(dx, -0.06, 0.045)
        path(c(-0.39, -0.16, 0.14, 0.37), c(0.43, 0.60, 0.60, 0.43),
             0.85, arrow = TRUE)
      } else if (hallmark == "Chronic inflammation") {
        a <- seq(0, 2 * pi, length.out = 19)
        radius <- rep(c(0.37, 0.22), length.out = length(a))
        polygon(-0.49 + radius * cos(a), 0.07 + radius * sin(a))
        dot(-0.49, 0.07, 0.10)
        ellipse(0.51, -0.12, 0.30)
        dot(0.51, -0.12, 0.11)
        for (q in list(c(-0.02, 0.28), c(0.13, 0.08), c(-0.05, -0.17))) {
          dot(q[1], q[2], 0.055)
        }
        path(c(-0.35, -0.08, 0.19, 0.42), c(0.58, 0.72, 0.66, 0.45),
             0.8, arrow = TRUE)
      } else if (hallmark == "Dysbiosis") {
        ellipse(-0.20, 0, 0.49, 0.29)
        for (dx in c(-0.46, -0.20, 0.06)) {
          path(c(dx, dx), c(0.30, 0.41), 0.65)
          path(c(dx, dx), c(-0.30, -0.41), 0.65)
        }
        t <- seq(0, 1, length.out = 70)
        path(0.29 + 0.58 * t, 0.15 * sin(2.4 * pi * t), 0.95)
        path(c(-0.47, -0.33, -0.20, -0.07, 0.06),
             c(0.03, 0.12, -0.10, 0.08, -0.04), 0.70)
      } else {
        stop("No icon is defined for hallmark: ", hallmark)
      }
      invisible(NULL)
    }
    
    parallel_sector_layout <- function(i, text_size) {
      info <- hallmark_protein_lines[[hallmarks12$hallmark[i]]]
      title_lines <- strsplit(hallmarks12$label[i], "\n", fixed = TRUE)[[1L]]
      if (hallmarks12$hallmark[i] == "Altered intercellular communication") {
        title_lines <- c("Altered intercellular", "communication")
      }
      active <- hallmarks12$hallmark[i] %in% represented_hallmarks
      if (active && length(title_lines) != 2L) {
        stop("A colored hallmark must have a two-line title: ", hallmarks12$hallmark[i])
      }
      labels <- if (active) {
        c(title_lines, paste(info$male, collapse = ", "), paste(info$female, collapse = ", "))
      } else {
        title_lines
      }
      colors <- if (active) c(COL$ink, COL$ink, COL$male, COL$female_text) else {
        rep(COL$Gray_text, length(labels))
      }
      n_lines <- length(labels)
      am <- (angs[i] + angs[i + 1L]) / 2
      rotation <- ((am * 180 / pi + 90) %% 180) - 90
      rotation <- round(rotation, 8)
      tr <- rotation * pi / 180
      dims <- vapply(labels, function(label) {
        measure_text(label, text_size, "bold", 1)
      }, numeric(2))
      em_mm <- text_size * 25.4 / 72
      padding_x <- 0.30
      padding_y <- 0.25
      ascent <- max(0.75 * em_mm, dims[2, ])
      descent <- 0.10 * em_mm
      line_pitch <- ascent + descent + 0.18
      block_width <- max(dims[1, ]) + 2 * padding_x
      block_height <- (n_lines - 1L) * line_pitch + ascent + descent + 2 * padding_y
      baseline_y <- padding_y + descent + rev(seq.int(0L, n_lines - 1L)) * line_pitch
      lo <- c(-block_width / 2, -block_height / 2)
      hi <- -lo
      corners <- expand.grid(x = c(lo[1], hi[1]), y = c(lo[2], hi[2]))
      for (radius in seq(from = 31.5, to = 40.5, by = 0.10)) {
        block_cx <- radius * cos(am)
        block_cy <- radius * sin(am)
        px <- block_cx + corners$x * cos(tr) - corners$y * sin(tr)
        py <- block_cy + corners$x * sin(tr) + corners$y * cos(tr)
        origin_x <- -(block_cx * cos(tr) + block_cy * sin(tr))
        origin_y <- block_cx * sin(tr) - block_cy * cos(tr)
        nearest <- c(max(lo[1], min(origin_x, hi[1])) - origin_x,
                     max(lo[2], min(origin_y, hi[2])) - origin_y)
        fits <- all(px^2 + py^2 <= (r1 - 0.7)^2) &&
          sum(nearest^2) >= text_inner_radius^2 &&
          all(cos(angs[i + 1L]) * py - sin(angs[i + 1L]) * px >= 0.7) &&
          all(sin(angs[i]) * px - cos(angs[i]) * py >= 0.7)
        if (fits) return(list(
          parallel_block = TRUE, labels = labels, colors = colors,
          fontsize = text_size, angle = rotation,
          cx = cx + block_cx, cy = cy + block_cy,
          width = block_width, height = block_height,
          left = padding_x, baseline_y = baseline_y, line_pitch = line_pitch
        ))
      }
      NULL
    }
    active_sectors <- which(hallmarks12$hallmark %in% represented_hallmarks)
    parallel_layouts <- vector("list", nrow(hallmarks12))
    common_font_size <- NA_real_
    for (candidate_size in seq(from = 7.2, to = 4.8, by = -0.1) * hallmark_scale) {
      trial <- lapply(active_sectors, function(i) parallel_sector_layout(i, candidate_size))
      if (all(vapply(trial, function(z) !is.null(z), logical(1)))) {
        parallel_layouts[active_sectors] <- trial
        common_font_size <- candidate_size
        break
      }
    }
    if (!is.finite(common_font_size)) {
      stop("The four-line hallmark blocks do not fit at the minimum readable font size.")
    }
    draw_parallel_sector <- function(z) {
      n_lines <- length(z$labels)
      if (!n_lines || any(grepl("\n", z$labels, fixed = TRUE)) ||
          length(z$fontsize) != 1L || !is.finite(z$fontsize) ||
          length(z$angle) != 1L || !is.finite(z$angle) ||
          length(z$baseline_y) != n_lines || length(z$colors) != n_lines ||
          any(abs(diff(z$baseline_y) + z$line_pitch) > 1e-10)) {
        stop("Invalid shared typography for the hallmark text block.")
      }
      grid::pushViewport(grid::viewport(
        x = u(z$cx), y = u(Y(z$cy)),
        width = grid::unit(z$width, "mm"), height = grid::unit(z$height, "mm"),
        angle = z$angle, clip = "off"
      ))
      on.exit(grid::popViewport(), add = TRUE)
      rows <- lapply(seq_len(n_lines), function(j) {
        grid::textGrob(
          label = z$labels[j],
          x = grid::unit(z$left, "mm"), y = grid::unit(z$baseline_y[j], "mm"),
          hjust = 0, vjust = 0, rot = 0, check.overlap = FALSE,
          gp = grid::gpar(
            fontfamily = resolved_font, fontsize = z$fontsize,
            fontface = "bold", col = z$colors[j], lineheight = 1
          )
        )
      })
      grid::grid.draw(grid::gTree(children = do.call(grid::gList, rows)))
      invisible(NULL)
    }
    sector_layout <- function(i) {
      if (hallmarks12$hallmark[i] %in% represented_hallmarks) {
        return(parallel_layouts[[i]])
      }
      for (factor in seq(from = 1, to = 0.58, length.out = 22L)) {
        z <- parallel_sector_layout(i, 7.4 * hallmark_scale * factor)
        if (!is.null(z)) return(z)
      }
      stop("Unable to fit text inside the hallmark sector: ", hallmarks12$hallmark[i])
    }
    layouts <- lapply(seq_len(nrow(hallmarks12)), sector_layout)
    for (i in seq_len(nrow(hallmarks12))) {
      a <- seq(angs[i], angs[i + 1L], length.out = 80)
      xx <- c(cx + r0 * cos(a), rev(cx + r1 * cos(a)))
      yy <- c(cy + r0 * sin(a), rev(cy + r1 * sin(a)))
      grid::grid.polygon(u(xx), u(Y(yy)), gp = grid::gpar(fill = hallmarks12$fill[i], col = "white", lwd = 1))
    }
    for (i in seq_len(nrow(hallmarks12))) {
      am <- (angs[i] + angs[i + 1L]) / 2
      icon_rotation <- round(((am * 180 / pi + 90) %% 180) - 90, 8)
      draw_hallmark_icon(
        hallmarks12$hallmark[i], cx + icon_radius * cos(am), cy + icon_radius * sin(am),
        size = icon_size, angle = icon_rotation,
        active = hallmarks12$hallmark[i] %in% represented_hallmarks
      )
    }
    circle(cx, cy, r0, "white", COL$rule, 0.8)
    text_at("Hallmarks\nof aging", cx, cy, 9.0 * hallmark_scale, COL$ink, "bold")
    for (z in layouts) draw_parallel_sector(z)
    
    # Re-draw inter-panel arrows on top so they are not hidden by panel content
    for (mid in separators) panel_arrow(mid, y = 80)
    
    invisible(NULL)
  }
  
  # Capture the unchanged, full-size drawing as Cairo SVG before shrinking it.
  # The SVG outlines its text. Thus fonts, line widths, symbols and arrows all
  # scale together, unlike calling draw() directly on a smaller graphics device.
  # The original four-panel draw() function and its large exports are unchanged.
  toc_picture <- NULL
  prepare_toc_picture <- function(svg_source = NULL) {
    if (!requireNamespace("grImport2", quietly = TRUE)) {
      stop("First run install.packages('grImport2').")
    }
    if (is.null(svg_source)) {
      svg_source <- tempfile(pattern = "nako_toc_source_", fileext = ".svg")
      on.exit(unlink(svg_source), add = TRUE)
      grDevices::svg(svg_source, width = 360 / 25.4, height = 150 / 25.4,
                     family = resolved_font, bg = "white", onefile = FALSE)
      device <- grDevices::dev.cur()
      tryCatch(draw(), finally = {
        devices <- grDevices::dev.list()
        if (!is.null(devices) && device %in% devices) grDevices::dev.off(device)
      })
    }
    if (!file.exists(svg_source)) stop("TOC source SVG not found: ", svg_source)
    grImport2::readPicture(svg_source, warn = FALSE)
  }

  draw_toc <- function() {
    if (is.null(toc_picture)) toc_picture <<- prepare_toc_picture()
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grImport2::grid.picture(
      toc_picture,
      x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
      width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
      just = "centre", distort = FALSE, expansion = 0, ext = "clipbbox"
    )
    invisible(NULL)
  }

  # Dedicated vector PDF export for the TOC; the main render_file() is unchanged.
  # All text is already outlined in the imported SVG, so no font substitution
  # occurs in the PDF. R rounds the PDF page dimensions down to whole points,
  # keeping the physical page within the requested 50 x 60 mm bounds.
  render_toc_file <- function(format) {
    basename_toc <- paste0(OUTPUT_BASENAME, "_TOC_50x60mm")
    if (format != "pdf") {
      return(render_file(format, draw_fun = draw_toc, basename_out = basename_toc,
                         width_mm = 50, height_mm = 60))
    }
    destination <- file.path(output_dir, paste0(basename_toc, ".pdf"))
    tmp <- tempfile(pattern = "nako_toc_", fileext = ".pdf")
    on.exit(unlink(tmp), add = TRUE)
    grDevices::pdf(tmp, width = 50 / 25.4, height = 60 / 25.4,
                   bg = "white", onefile = FALSE, useDingbats = FALSE,
                   compress = TRUE, version = "1.4", title = "")
    device <- grDevices::dev.cur()
    tryCatch(draw_toc(), finally = {
      devices <- grDevices::dev.list()
      if (!is.null(devices) && device %in% devices) grDevices::dev.off(device)
    })
    if (!file.exists(tmp) || !is.finite(file.info(tmp)$size) || file.info(tmp)$size <= 0) {
      stop("TOC PDF export failed.")
    }
    if (!file.copy(tmp, destination, overwrite = TRUE) ||
        file.info(destination)$size != file.info(tmp)$size) stop("Unable to save: ", destination)
    if (file.info(destination)$size > 10 * 1024^2) warning("Export exceeds 10 MB: ", basename(destination), call. = FALSE)
    message("Saved: ", destination)
    destination
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  render_file <- function(format, draw_fun = draw, basename_out = OUTPUT_BASENAME,
                          width_mm = 360, height_mm = 150) {
    destination <- file.path(output_dir, paste0(basename_out, ".", format))
    tmp <- tempfile(pattern = "nako_graphical_", fileext = paste0(".", format))
    on.exit(unlink(tmp), add = TRUE)
    if (format == "pdf") {
      grDevices::cairo_pdf(tmp, width = width_mm / 25.4, height = height_mm / 25.4,
                           family = resolved_font, bg = "white", onefile = FALSE, fallback_resolution = DPI)
    } else if (format == "svg") {
      grDevices::svg(tmp, width = width_mm / 25.4, height = height_mm / 25.4,
                     family = resolved_font, bg = "white", onefile = FALSE)
    } else {
      f <- switch(format, png = grDevices::png, jpg = grDevices::jpeg, tiff = grDevices::tiff)
      if (is.null(f)) stop("Unsupported export format: ", format)
      args <- list(filename = tmp, width = width_mm, height = height_mm, units = "mm", res = DPI,
                   type = "cairo", family = resolved_font, bg = "white", antialias = "gray")
      if (format == "jpg") args$quality <- 100
      if (format == "tiff") args$compression <- "lzw"
      do.call(f, args)
    }
    device <- grDevices::dev.cur()
    tryCatch(draw_fun(), finally = {
      devices <- grDevices::dev.list()
      if (!is.null(devices) && device %in% devices) grDevices::dev.off(device)
    })
    if (!file.exists(tmp) || !is.finite(file.info(tmp)$size) || file.info(tmp)$size <= 0) {
      stop("Graphics export failed: ", format)
    }
    if (!file.copy(tmp, destination, overwrite = TRUE) ||
        file.info(destination)$size != file.info(tmp)$size) stop("Unable to save: ", destination)
    if (file.info(destination)$size > 10 * 1024^2) warning("Export exceeds 10 MB: ", basename(destination), call. = FALSE)
    message("Saved: ", destination)
    destination
  }
  
  formats <- c("pdf", "png", "tiff", "jpg")
  if (EXPORT_SVG) formats <- c(formats, "svg")
  written <- vapply(formats, render_file, character(1))
  toc_files <- character()
  if (EXPORT_TOC) {
    # Reuse the freshly exported full-size SVG; otherwise create a temporary one.
    # Both routes use the unchanged draw() at its original 360 x 150 mm size.
    toc_picture <- prepare_toc_picture(if (EXPORT_SVG) unname(written["svg"]) else NULL)
    toc_files <- vapply(c("pdf", "png", "tiff"), render_toc_file, character(1))
    message("TOC: complete figure, approximately 50 x 20.83 mm, centered on a 50 x 60 mm page.")
    warning("The uniformly reduced TOC retains all content, but its text is very small. Check legibility at actual size before submission.", call. = FALSE)
  }
  
  report <- data.frame(
    item = c("Initially sampled", "Initial male", "Initial female", "Excluded after QC",
             "Participants", "Male", "Female", "Proteins quantified", "Proteins retained",
             "Male signature", "Female signature", "Shared proteins", "Shared symbols",
             "AgeProt R2 male", "AgeProt R2 female", "Font"),
    value = c(initial_participants, initial_n_male, initial_n_female, excluded_after_qc,
              ncol(expr), n_male, n_female, INITIAL_PROTEINS, nrow(expr),
              length(sets$male), length(sets$female), length(shared),
              paste(symbols(shared), collapse = "; "), sprintf("%.6f", ageprot_r2["male"]),
              sprintf("%.6f", ageprot_r2["female"]), resolved_font)
  )
  utils::write.csv(report, file.path(output_dir, "graphical_abstract_summary.csv"), row.names = FALSE)
  utils::write.csv(as.data.frame.matrix(quota), file.path(output_dir, "graphical_abstract_quota_counts.csv"))
  utils::write.csv(as.data.frame(initial_quota),
                   file.path(output_dir, "graphical_abstract_initial_quota_counts.csv"))
  utils::write.csv(data.frame(hallmark = names(hallmark_protein_lines), male = vapply(hallmark_protein_lines, function(x) paste(x$male, collapse = "; "), character(1)), female = vapply(hallmark_protein_lines, function(x) paste(x$female, collapse = "; "), character(1))), file.path(output_dir, "graphical_abstract_hallmark_proteins.csv"), row.names = FALSE)
  manifest <- data.frame(file = basename(inputs), md5 = unname(tools::md5sum(inputs)))
  utils::write.csv(manifest, file.path(output_dir, "graphical_abstract_inputs.csv"), row.names = FALSE)
  capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo_graphical_abstract.txt"))
  if (SHOW_IN_PLOT_PANE) draw()
  toc_summary <- sprintf(
    paste0("Of %d sampled NAKO-Bremen participants, %d passed QC. Signatures comprised %d male proteins (%s) and %d female proteins (%s). ",
           "%s were shared. ",
           "Chronological age explained %.0f%% and %.0f%% of AgeProt variance in male and female participants, respectively."),
    initial_participants, ncol(expr), length(sets$male), paste(sort(selected_symbols$male), collapse = ", "),
    length(sets$female), paste(sort(selected_symbols$female), collapse = ", "),
    paste(symbols(shared), collapse = " and "),
    100 * ageprot_r2["male"], 100 * ageprot_r2["female"]
  )
  toc_word_count <- length(strsplit(trimws(toc_summary), "[[:space:]]+")[[1L]])
  if (toc_word_count > 80L) warning("The TOC text exceeds 80 words: ", toc_word_count, call. = FALSE)
  writeLines(c(TITLE, "", toc_summary), file.path(output_dir, "graphical_abstract_TOC_text_draft.txt"), useBytes = TRUE)
  figure_legend <- c(
    paste0("Graphical abstract. ", TITLE, "."),
    sprintf(
      paste0("(1) Study design: initial quota sampling of %d NAKO-Bremen participants (%d male and %d female), ",
             "with 12 participants in each age-by-sex stratum. ",
             "(2) Analytical pipeline: LC-MS/MS profiling and post-LC-MS/MS quality control, ",
             "with %d participants retained; imputation and normalization yielded %d proteins for analysis, ",
             "followed by sex-stratified ordinal elastic net regression, ",
             "unsupervised internal assessment, and construction of the AgeProt score. ",
             "(3) Age-associated signatures: volcano plots show limma age coefficients per 10 years and nominal p-values; ",
             "the Venn diagram shows the overlap between the selected signatures. ",
             "Colored boxes list all selected proteins for each sex, including proteins shared between the signatures. ",
             "Chronological age explained %.0f%% and %.0f%% of AgeProt variance in male and female participants, respectively. ",
             "(4) Biological interpretation: descriptive mapping to the Hallmarks of Aging. ",
             "Green and orange identify male and female participants throughout; the five represented hallmarks occupy consecutive blue sectors, whereas gray sectors are not represented in the current mapping. The small vector icons are schematic illustrations of each hallmark."),
      initial_participants, initial_n_male, initial_n_female, ncol(expr), nrow(expr),
      100 * ageprot_r2["male"], 100 * ageprot_r2["female"]
    ),
    sprintf("Selected proteins in male participants (n=%d): %s.", length(sets$male),
            paste(sort(selected_symbols$male), collapse = ", ")),
    sprintf("Selected proteins in female participants (n=%d): %s.", length(sets$female),
            paste(sort(selected_symbols$female), collapse = ", ")),
    paste0("Shared proteins: ", paste(symbols(shared), collapse = " and "), ".")
  )
  legend_file <- file.path(output_dir, "graphical_abstract_figure_legend.txt")
  writeLines(figure_legend, legend_file, useBytes = TRUE)
  writeLines(c(
    "Aging Cell graphical abstract: export notes",
    "Updated layout: compact panels 1-2, centered headings, identical inter-panel arrows, enlarged hallmark wheel.",
    "Panel 1 shows the initial quota sample; panel 2 shows the analytical sample retained after QC.",
    "The initial quotas are design counts and do not add observations to the analytical data.",
    "The four panels use the same bottom alignment for their final text blocks.",
    "Venn boxes contain every selected protein once, arranged in two columns for each sex.",
    "Volcano panels have matching physical dimensions and retain their original sex-specific x- and y-axis scales.",
    "Hallmark protein labels have no background patches; text retains its shared rotation and baselines.",
    "Hallmark icons are drawn directly with R grid primitives; no external images or icon packages are required.",
    "All five represented hallmarks have the same blue fill and occupy adjacent sectors.",
    "The icon vocabulary was informed by Figure 1 of Lopez-Otin et al. (2023), Cell, doi:10.1016/j.cell.2022.11.001.",
    paste("Font:", resolved_font),
    paste("Raster resolution:", DPI, "dpi; Cairo graphics devices."),
    "The 360 x 150 mm output is the detailed central illustration, not a size-compliant TOC image.",
    "The separate 50 x 60 mm TOC page contains the complete four-panel illustration, uniformly scaled and centered.",
    "The illustration occupies approximately 50 x 20.83 mm; the remaining space is white to avoid stretching or cropping.",
    "TOC conversion uses grImport2 and Cairo SVG outlines; the PDF retains vector paths, while PNG and TIFF use 600 dpi.",
    "No panels or labels are removed. Small text may be illegible at actual size; size compliance alone does not ensure legibility.",
    "The full illustration retains all twelve hallmarks. Gray means not represented in the current descriptive mapping.",
    "Statistical models, age values and saved protein selections have not been redefined by this graphics script.",
    paste("TOC summary word count:", toc_word_count),
    "Both the TOC text and the separate full figure legend include every selected protein symbol.",
    "The TOC text draft must be accompanied by the final author list with the corresponding author marked by an asterisk.",
    "Check the final title and author list against the submission manuscript.",
    "Author guidelines, accessed 28 September 2026:",
    "https://onlinelibrary.wiley.com/page/journal/14749726/homepage/forauthors.html",
    "Wiley figure guidance:",
    "https://authors.wiley.com/author-resources/Journal-Authors/Prepare/manuscript-preparation-guidelines.html/figure-preparation.html"
  ), file.path(output_dir, "graphical_abstract_export_notes.txt"), useBytes = TRUE)
  message("\nCompleted. Output: ", output_dir)
  invisible(list(output_dir = output_dir, files = written, toc_files = toc_files, summary = report,
                 quota = quota, initial_quota = initial_quota, ageprot_r2 = ageprot_r2,
                 hallmark_proteins = hallmark_protein_lines,
                 venn_box_proteins = venn_box_proteins, legend_file = legend_file, toc_word_count = toc_word_count,
                 draw = draw, draw_toc = draw_toc))
})
