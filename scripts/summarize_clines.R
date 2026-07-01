########################################################################################
# TREND GROUPING, QC FILTERING, FINAL PLOTS, AND PAPER TABLES
########################################################################################


# Usage example:
#   Rscript summarize_clines_single.R
#
# Expects fitted-results table from the fitting script:
#   hzar_ML_summary_with_CI.tsv

out_hzar_dir <- "hzar_inputs_all"

# Single input table from cline_fitting_hzar_parboot_single.R
in_param_tsv <- "hzar_ML_summary.tsv"

# Clean output directory. Change to "." if you prefer outputs in the working directory.
out_dir <- "cline_results"
fs::dir_create(out_dir)

# Main reproducible outputs.
# 1) all fitted loci before QC, including trend, metadata, and QC flags
# 2) the same table after QC filtering
out_results_pre_qc_tsv <- fs::path(out_dir, "hzar_results_pre_QC_full.tsv")
out_results_qc_tsv     <- fs::path(out_dir, "hzar_results_QC_filtered_full.tsv")

# Paper-facing outputs.
# Table 2-style group summary: medians with IQRs in brackets.
# Table S6-style per-locus table: point estimates with 95% parametric-bootstrap CIs in brackets.
out_group_paper_tsv <- fs::path(out_dir, "cline_summary_by_trend.tsv")
out_locus_paper_tsv <- fs::path(out_dir, "cline_per_locus_parameters_QC_filtered.tsv")

# Plot output.
out_violin_pdf <- fs::path(out_dir, "Hzar_center_width_violin.pdf")

# Small text manifest explaining the output files.
out_manifest_txt <- fs::path(out_dir, "README_outputs.txt")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(patchwork)
})

write_tsv_plain <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

# Numeric output columns expected from the HZAR fitting step.
# Includes columns from both older site-bootstrap and newer parametric-bootstrap outputs.
numeric_fit_cols <- c(
  "n_sites", "alt_min", "alt_max", "alt_rng", "dp_obs",
  "center", "center_lo", "center_hi",
  "width", "width_lo", "width_hi",
  "pMin", "pMax", "p_low_orig", "p_high_orig",
  "delta_p", "delta_p_signed", "delta_p_abs",
  "delta_p_abs_lo", "delta_p_abs_hi",
  "logLik", "AIC",
  "rho", "rho_spearman",
  "center_bound_low", "center_bound_high", "width_bound_high",
  "n_boot_total", "n_boot_ok", "boot_success_prop",
  "center_boot_sd", "width_boot_sd", "delta_p_abs_boot_sd"
)

coerce_fit_numeric_cols <- function(x) {
  x %>% mutate(across(any_of(numeric_fit_cols), as.numeric))
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

ensure_numeric_column <- function(x, col) {
  if (!col %in% names(x)) x[[col]] <- NA_real_
  x[[col]] <- suppressWarnings(as.numeric(x[[col]]))
  x
}

ensure_character_column <- function(x, col) {
  if (!col %in% names(x)) x[[col]] <- NA_character_
  x[[col]] <- as.character(x[[col]])
  x
}

format_median_iqr <- function(med, q1, q3, digits = 0) {
  fmt <- function(z) formatC(z, format = "f", digits = digits)
  ifelse(
    is.finite(med) & is.finite(q1) & is.finite(q3),
    paste0(fmt(med), " [", fmt(q1), "–", fmt(q3), "]"),
    "NA"
  )
}

format_est_ci <- function(est, lo, hi, digits = 2) {
  fmt <- function(z) formatC(z, format = "f", digits = digits)
  ifelse(
    is.finite(est) & is.finite(lo) & is.finite(hi),
    paste0(fmt(est), " [", fmt(lo), "–", fmt(hi), "]"),
    "NA"
  )
}

trend_from_file <- function(path_tsv, rho_threshold = 0.10) {
  df <- tryCatch(
    read.table(path_tsv, header = TRUE, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)

  df <- df %>%
    mutate(
      altitude = as.numeric(altitude),
      k_minor  = as.numeric(k_minor),
      n_chrom  = as.numeric(n_chrom),
      p        = k_minor / n_chrom
    ) %>%
    filter(is.finite(altitude), is.finite(p), is.finite(n_chrom), n_chrom > 0)

  if (nrow(df) < 3) return(tibble(trend = "Flat", rho = NA_real_))

  rho <- suppressWarnings(cor(df$altitude, df$p, method = "spearman"))
  trend <- case_when(
    is.finite(rho) & rho >  rho_threshold ~ "Increasing",
    is.finite(rho) & rho < -rho_threshold ~ "Decreasing",
    TRUE                                  ~ "Flat"
  )

  tibble(trend = trend, rho = rho)
}

## 1. Read the single HZAR fitted-results table
if (file_exists(in_param_tsv)) {
  message("Reading fitted HZAR results from: ", in_param_tsv)
  res <- read_tsv(in_param_tsv, col_types = cols(.default = col_character())) %>%
    distinct(locus_id, .keep_all = TRUE) %>%
    coerce_fit_numeric_cols()
} else {
  stop("Error: fitted-results table not found: ", in_param_tsv)
}


if ("boot_success_prop" %in% names(res) && any(is.finite(res$boot_success_prop))) {
  boot_summary <- res %>%
    summarise(
      n_loci_with_bootstrap_info = sum(is.finite(boot_success_prop)),
      median_boot_success_prop = median(boot_success_prop, na.rm = TRUE),
      min_boot_success_prop = min(boot_success_prop, na.rm = TRUE),
      n_loci_boot_success_lt_0_8 = sum(is.finite(boot_success_prop) & boot_success_prop < 0.80)
    )
  message(
    "Bootstrap diagnostics: median success proportion = ",
    round(boot_summary$median_boot_success_prop, 3),
    "; minimum = ", round(boot_summary$min_boot_success_prop, 3),
    "; loci < 0.80 = ", boot_summary$n_loci_boot_success_lt_0_8
  )
}

if ("bootstrap_type" %in% names(res)) {
  message("Bootstrap type(s): ", paste(unique(res$bootstrap_type), collapse = ", "))
}
if (all(c("center_ci_hits_bound", "width_ci_hits_bound") %in% names(res))) {
  message(
    "CI bound diagnostics: center bound hits = ",
    sum(as_logical_safe(res$center_ci_hits_bound), na.rm = TRUE),
    "; width bound hits = ",
    sum(as_logical_safe(res$width_ci_hits_bound), na.rm = TRUE)
  )
}

## 2. Calculate Spearman rho trends from per-locus HZAR input files when available
paths <- fs::path(out_hzar_dir, paste0(gsub("[^A-Za-z0-9_.:-]+", "_", res$locus_id), ".tsv"))
names(paths) <- res$locus_id

trend_tbl <- bind_rows(lapply(names(paths), function(lid) {
  out <- trend_from_file(paths[[lid]])
  if (!is.null(out)) mutate(out, locus_id = lid)
}))

# Preserve fitted-script amplitude columns if they are already present.
# Older fitting outputs do not contain delta_p_abs, so compute it here as a fallback.
if (!"delta_p" %in% names(res)) {
  res <- res %>% mutate(delta_p = pMax - pMin)
}
if (!"delta_p_abs" %in% names(res)) {
  res <- res %>% mutate(delta_p_abs = abs(delta_p))
}

# Avoid trend/rho suffixes if re-running this script on previously merged output.
res <- res %>% select(-any_of(c("trend", "rho")))

if (nrow(trend_tbl) > 0) {
  res_join <- res %>%
    left_join(trend_tbl, by = "locus_id")
} else {
  warning("No per-locus HZAR TSVs found for trend recalculation; falling back to trend_spearman/rho_spearman when available.")
  res_join <- res
}

# Fallback to trend/rho already emitted by the fitting script.
res_join <- res_join %>%
  ensure_character_column("trend") %>%
  ensure_numeric_column("rho") %>%
  ensure_character_column("trend_spearman") %>%
  ensure_numeric_column("rho_spearman") %>%
  mutate(
    trend = coalesce(trend, trend_spearman),
    rho   = coalesce(rho, rho_spearman)
  ) %>%
  arrange(locus_id)


## 3. QC filtering for unreliable loci

QC_RANGE_BUFFER_PROP <- 0.10     # allow center within ±10% of sampled range
QC_WIDTH_MAX_FACTOR  <- 5        # max width allowed as multiple of sampled range
QC_WIDTH_MIN_ABS     <- 10       # widths < 10 m are too sharp for site spacing
QC_DP_MIN            <- 0.02     # minimal fitted amplitude
QC_DP_OBS_MIN        <- 0.02     # minimal observed amplitude

get_locus_meta <- function(locus_id) {
  f <- fs::path(out_hzar_dir, paste0(gsub("[^A-Za-z0-9_.:-]+", "_", locus_id), ".tsv"))
  if (!file.exists(f)) return(NULL)

  df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)

  df <- df %>%
    mutate(
      altitude = as.numeric(altitude),
      k_minor  = as.numeric(k_minor),
      n_chrom  = as.numeric(n_chrom),
      p        = k_minor / n_chrom
    ) %>%
    filter(is.finite(altitude), is.finite(p), is.finite(n_chrom), n_chrom > 0)

  if (nrow(df) < 3) return(NULL)

  tibble(
    locus_id = locus_id,
    alt_min  = min(df$altitude, na.rm = TRUE),
    alt_max  = max(df$altitude, na.rm = TRUE),
    alt_rng  = max(df$altitude, na.rm = TRUE) - min(df$altitude, na.rm = TRUE),
    dp_obs   = diff(range(df$p, na.rm = TRUE))
  )
}

# Per-locus metadata may already be present when using the parametric-bootstrap fitting script.
# File-derived metadata is used only to fill missing values.
meta_cols <- c("alt_min", "alt_max", "alt_rng", "dp_obs")

meta_list <- lapply(res_join$locus_id, get_locus_meta)
meta_list <- meta_list[!vapply(meta_list, is.null, logical(1))]

meta_tbl <- if (length(meta_list) > 0) {
  bind_rows(meta_list) %>%
    mutate(across(all_of(meta_cols), as.numeric))
} else {
  tibble(
    locus_id = character(),
    alt_min  = numeric(),
    alt_max  = numeric(),
    alt_rng  = numeric(),
    dp_obs   = numeric()
  )
}

res_qc_input <- res_join %>%
  left_join(
    meta_tbl %>% rename_with(~ paste0(.x, "_from_file"), all_of(meta_cols)),
    by = "locus_id"
  )

for (col in meta_cols) {
  file_col <- paste0(col, "_from_file")
  res_qc_input <- ensure_numeric_column(res_qc_input, col)
  if (file_col %in% names(res_qc_input)) {
    res_qc_input[[file_col]] <- suppressWarnings(as.numeric(res_qc_input[[file_col]]))
    res_qc_input[[col]] <- dplyr::coalesce(res_qc_input[[col]], res_qc_input[[file_col]])
  }
}

res_qc_input <- res_qc_input %>% select(-any_of(paste0(meta_cols, "_from_file")))

missing_meta_n <- res_qc_input %>%
  summarise(n_missing_meta = sum(!is.finite(alt_min) | !is.finite(alt_max) | !is.finite(alt_rng) | !is.finite(dp_obs))) %>%
  pull(n_missing_meta)

if (missing_meta_n > 0) {
  warning(missing_meta_n, " loci lack complete altitude metadata and will be flagged by QC.")
}

res_qc <- res_qc_input %>%
  mutate(
    flag_center_extrap = ifelse(
      is.finite(center) & is.finite(alt_min) & is.finite(alt_max) & is.finite(alt_rng),
      center < (alt_min - QC_RANGE_BUFFER_PROP * alt_rng) |
        center > (alt_max + QC_RANGE_BUFFER_PROP * alt_rng),
      TRUE
    ),
    flag_width_inflate = ifelse(
      is.finite(width) & is.finite(alt_rng),
      width > QC_WIDTH_MAX_FACTOR * alt_rng | width < QC_WIDTH_MIN_ABS,
      TRUE
    ),
    flag_delta_p_small = !is.finite(delta_p_abs) | delta_p_abs < QC_DP_MIN,
    flag_dp_obs_small  = !is.finite(dp_obs)  | dp_obs  < QC_DP_OBS_MIN,
    flagged = flag_center_extrap | flag_width_inflate | flag_delta_p_small | flag_dp_obs_small
  )

write_tsv_plain(res_qc %>% arrange(locus_id), out_results_pre_qc_tsv)
message("Wrote: ", out_results_pre_qc_tsv)

res_qc_filt <- res_qc %>%
  filter(!flagged, is.finite(center), is.finite(width), width > 0) %>%
  arrange(locus_id)

write_tsv_plain(res_qc_filt, out_results_qc_tsv)
message("Wrote: ", out_results_qc_tsv)

message(
  "QC kept ", nrow(res_qc_filt), " / ", nrow(res_qc), " loci (",
  round(100 * nrow(res_qc_filt) / max(1, nrow(res_qc)), 1), "%)."
)

## 4. Generate box/violin visualization metrics

res_plot <- res_qc_filt %>%
  filter(is.finite(center), is.finite(width), width > 0, trend %in% c("Increasing", "Decreasing", "Flat"))

p_center_box <- ggplot(res_plot, aes(trend, center, fill = trend)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c(Increasing = "#1b9e77", Decreasing = "#d95f02", Flat = "grey60")) +
  labs(x = NULL, y = "Cline center (m, log10 scale)", fill = "Trend") +
  theme_minimal(base_size = 11)

p_width_box <- ggplot(res_plot, aes(trend, width, fill = trend)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c(Increasing = "#1b9e77", Decreasing = "#d95f02", Flat = "grey60")) +
  labs(x = NULL, y = "Cline width (m, log10 scale)", fill = "Trend") +
  theme_minimal(base_size = 11)

ggsave(out_violin_pdf, (p_center_box | p_width_box), width = 180, height = 95, units = "mm")
message("Wrote vector plot layout to: ", out_violin_pdf)

## 5. Build paper-facing and detailed group summaries

summary_dat <- res_qc_filt %>%
  filter(trend %in% c("Increasing", "Decreasing")) %>%
  mutate(
    trend = factor(trend, levels = c("Decreasing", "Increasing")),
    center_ci_hits_bound_logical = if ("center_ci_hits_bound" %in% names(.)) as_logical_safe(center_ci_hits_bound) else FALSE,
    width_ci_hits_bound_logical  = if ("width_ci_hits_bound"  %in% names(.)) as_logical_safe(width_ci_hits_bound)  else FALSE
  )

group_detailed <- summary_dat %>%
  group_by(trend) %>%
  summarise(
    n_loci = n(),
    center_med = median(center, na.rm = TRUE),
    center_q1  = quantile(center, 0.25, na.rm = TRUE),
    center_q3  = quantile(center, 0.75, na.rm = TRUE),
    center_IQR = IQR(center, na.rm = TRUE),
    width_med  = median(width, na.rm = TRUE),
    width_q1   = quantile(width, 0.25, na.rm = TRUE),
    width_q3   = quantile(width, 0.75, na.rm = TRUE),
    width_IQR  = IQR(width, na.rm = TRUE),
    dP_abs_med = median(delta_p_abs, na.rm = TRUE),
    dP_abs_q1  = quantile(delta_p_abs, 0.25, na.rm = TRUE),
    dP_abs_q3  = quantile(delta_p_abs, 0.75, na.rm = TRUE),
    dP_abs_IQR = IQR(delta_p_abs, na.rm = TRUE),
    rho_med    = median(rho, na.rm = TRUE),
    rho_q1     = quantile(rho, 0.25, na.rm = TRUE),
    rho_q3     = quantile(rho, 0.75, na.rm = TRUE),
    rho_IQR    = IQR(rho, na.rm = TRUE),
    AIC_med    = median(AIC, na.rm = TRUE),
    boot_success_prop_med = if ("boot_success_prop" %in% names(summary_dat)) median(boot_success_prop, na.rm = TRUE) else NA_real_,
    center_ci_hits_bound_n = sum(center_ci_hits_bound_logical, na.rm = TRUE),
    width_ci_hits_bound_n  = sum(width_ci_hits_bound_logical, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(trend)

# Paper-ready Table 2 style summary: medians with Q1-Q3 in brackets.
group_paper <- group_detailed %>%
  transmute(
    Trend = as.character(trend),
    N_loci = n_loci,
    Center_m = format_median_iqr(center_med, center_q1, center_q3, digits = 0),
    Width_m = format_median_iqr(width_med, width_q1, width_q3, digits = 0),
    Amplitude_abs_delta_p = format_median_iqr(dP_abs_med, dP_abs_q1, dP_abs_q3, digits = 2),
    Spearman_rho = format_median_iqr(rho_med, rho_q1, rho_q3, digits = 2)
  )

write_tsv_plain(group_paper, out_group_paper_tsv)
message("Wrote: ", out_group_paper_tsv)

## 6. Build paper-facing per-locus supplementary table

# Use original-frequency columns when available. These are preferable to pMin/pMax for reader-facing output
# because direction-aware fitting may flip decreasing loci internally.
# The reader-facing table reports point estimates with 95% parametric-bootstrap CIs in brackets.
# A separate components table is also written for machine-readable access to the individual CI bounds.
locus_paper_components <- res_qc_filt %>%
  ensure_numeric_column("rho") %>%
  ensure_numeric_column("n_sites") %>%
  ensure_numeric_column("center") %>%
  ensure_numeric_column("center_lo") %>%
  ensure_numeric_column("center_hi") %>%
  ensure_numeric_column("width") %>%
  ensure_numeric_column("width_lo") %>%
  ensure_numeric_column("width_hi") %>%
  ensure_numeric_column("p_low_orig") %>%
  ensure_numeric_column("p_high_orig") %>%
  ensure_numeric_column("delta_p_signed") %>%
  ensure_numeric_column("delta_p_abs") %>%
  ensure_numeric_column("delta_p_abs_lo") %>%
  ensure_numeric_column("delta_p_abs_hi") %>%
  mutate(
    # fallback for older outputs without original-orientation columns
    p_low_orig       = coalesce(p_low_orig, pMin),
    p_high_orig      = coalesce(p_high_orig, pMax),
    delta_p_signed   = coalesce(delta_p_signed, delta_p),
    delta_p_abs      = coalesce(delta_p_abs, abs(delta_p_signed)),
    delta_p_abs_lo   = coalesce(delta_p_abs_lo, NA_real_),
    delta_p_abs_hi   = coalesce(delta_p_abs_hi, NA_real_),
    trend = factor(trend, levels = c("Decreasing", "Increasing", "Flat"))
  ) %>%
  arrange(trend, locus_id) %>%
  transmute(
    Locus_ID = locus_id,
    Trend = as.character(trend),
    N_sites = as.integer(round(n_sites)),
    Spearman_rho = round(rho, 2),
    Center = round(center, 2),
    Center_CI_low = round(center_lo, 2),
    Center_CI_high = round(center_hi, 2),
    Width = round(width, 2),
    Width_CI_low = round(width_lo, 2),
    Width_CI_high = round(width_hi, 2),
    p_low_altitude = round(p_low_orig, 3),
    p_high_altitude = round(p_high_orig, 3),
    Amplitude_abs_delta_p = round(delta_p_abs, 3),
    Amplitude_CI_low = round(delta_p_abs_lo, 3),
    Amplitude_CI_high = round(delta_p_abs_hi, 3)
  )

locus_paper <- locus_paper_components %>%
  transmute(
    Locus_ID,
    Trend,
    N_sites,
    Spearman_rho,
    Center_m_95CI = format_est_ci(Center, Center_CI_low, Center_CI_high, digits = 2),
    Width_m_95CI = format_est_ci(Width, Width_CI_low, Width_CI_high, digits = 2),
    p_low_altitude,
    p_high_altitude,
    Amplitude_abs_delta_p_95CI = format_est_ci(
      Amplitude_abs_delta_p,
      Amplitude_CI_low,
      Amplitude_CI_high,
      digits = 3
    )
  )

write_tsv_plain(locus_paper, out_locus_paper_tsv)
message("Wrote: ", out_locus_paper_tsv)


writeLines(
  c(
    "Cline analysis output files",
    "===========================",
    "",
    "hzar_results_pre_QC_full.tsv",
    "  Full merged HZAR results before QC. Includes fitted parameters, bootstrap diagnostics, trend columns, altitude metadata, and QC flags.",
    "",
    "hzar_results_QC_filtered_full.tsv",
    "  Full merged HZAR results after QC filtering. This is the main machine-readable results table for downstream checks.",
    "",
    "cline_summary_by_trend.tsv",
    "  Main summary table. Values are medians with interquartile ranges in brackets.",
    "",
    "cline_per_locus_parameters_QC_filtered_paper.tsv",
    "  Main per-locus table. Center, width, and amplitude are formatted as point estimate [95% parametric-bootstrap CI].",
    "",
    "Hzar_center_width_violin.pdf",
    "  Distribution of fitted cline centers and widths by trend group, using QC-filtered loci.",
    ""
  ),
  con = out_manifest_txt
)
message("Wrote: ", out_manifest_txt)

print(table(res_plot$trend))
