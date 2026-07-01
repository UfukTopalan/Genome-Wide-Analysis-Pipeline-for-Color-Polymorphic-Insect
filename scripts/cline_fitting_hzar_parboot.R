# Authors: Ufuk Topalan, Ismail Kudret Saglam
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------

########################################################################################
# hzar input generation and single-output cline fitting
# direction-aware fitting + parametric binomial bootstrap
########################################################################################


# Usage examples:
#   Rscript cline_fitting_hzar_parboot_single.R --cores=4 --n_boot=100 --force_rebuild_inputs
#   Rscript cline_fitting_hzar_parboot_single.R --cores=2 --no_bootstrap --max_loci=5 --force_rebuild_inputs
#
# This simple version does not use Slurm arrays and does not split loci into chunks.
# It writes one fitted-results table: hzar_ML_summary_with_CI.tsv

## ====== configuration ================================================================
geno_path              <- "isophya71.assoc.geno.tsv"
info_path              <- "isophya71.info"
recode_to_global_minor <- TRUE
min_sites_per_locus    <- 10L
out_hzar_dir           <- "hzar_inputs_all"
out_param_tsv          <- "hzar_ML_summary.tsv"
out_input_log_tsv      <- "hzar_input_generation_log.tsv"

# hzar model options
do_bootstrap           <- TRUE
bootstrap_type         <- "parametric"  # "parametric" or "site"
n_boot                 <- 100L

# interpretation controls
orient_to_increasing   <- TRUE           # decreasing loci are fitted as 1 - p for comparable center/width
rho_threshold          <- 0.10

# parameter-space bounds used during fitting
# these match the downstream qc philosophy and prevent the optimizer from exploring uninformative regions.
use_parameter_bounds   <- TRUE
center_buffer_prop     <- 0.10           # allow center within sampled range +/- 10% of range
width_max_factor       <- 5              # maximum width as multiple of sampled altitude range

# main fit and bootstrap fit controls
main_burnin            <- 10000L
main_mcmc              <- 1000000L
boot_burnin            <- 5000L
boot_mcmc              <- 50000L
base_seed              <- 1L

## ====== command-line arguments =======================================================
args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}

has_flag <- function(flag) {
  flag %in% args
}

core_arg <- get_arg_value("--cores")
if (!is.null(core_arg)) {
  n_cores <- as.integer(core_arg)
} else {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
  n_cores <- if (slurm_cpus != "") as.integer(slurm_cpus) else max(1L, parallel::detectCores(logical = TRUE) - 1L)
}
n_cores <- max(1L, n_cores)

n_boot_arg <- get_arg_value("--n_boot")
if (!is.null(n_boot_arg)) n_boot <- as.integer(n_boot_arg)


bootstrap_type_arg <- get_arg_value("--bootstrap_type")
if (!is.null(bootstrap_type_arg)) bootstrap_type <- bootstrap_type_arg
if (!bootstrap_type %in% c("parametric", "site")) stop("--bootstrap_type must be 'parametric' or 'site'")

if (has_flag("--no_bootstrap")) do_bootstrap <- FALSE
if (has_flag("--no_orientation")) orient_to_increasing <- FALSE
if (has_flag("--no_bounds")) use_parameter_bounds <- FALSE

force_rebuild_inputs <- has_flag("--force_rebuild_inputs")
max_loci_arg <- get_arg_value("--max_loci")
max_loci <- if (is.null(max_loci_arg)) NA_integer_ else as.integer(max_loci_arg)

## ====== libraries ====================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(fs)
  library(hzar)
  library(parallel)
})

## ====== prevent thread over-subscription ============================================
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
}

## ====== reproducible rng ============================================================
RNGkind("L'Ecuyer-CMRG")
set.seed(base_seed)

message("Using ", n_cores, " CPU core(s) for parallel processing.")
message("Bootstrap enabled: ", do_bootstrap, " | type = ", bootstrap_type, " | n_boot = ", n_boot)
message("Direction-aware fitting: ", orient_to_increasing)
message("Parameter bounds during fitting: ", use_parameter_bounds)
message("Seed = ", base_seed)

## ====== helpers ======================================================================
convert_to_minor_counts <- function(gt_chr, maj, min) {
  gt_chr <- toupper(trimws(gt_chr))
  out <- rep(NA_integer_, length(gt_chr))
  mm <- paste0(maj, maj); mn <- paste0(maj, min)
  nm <- paste0(min, maj); nn <- paste0(min, min)
  out[gt_chr == mm] <- 0L
  out[gt_chr == mn] <- 1L
  out[gt_chr == nm] <- 1L
  out[gt_chr == nn] <- 2L
  out[gt_chr %in% c("NN", "NA", ".", "")] <- NA_integer_
  out
}

write_tsv_plain <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

sanitize_locus_id <- function(x) {
  gsub("[^A-Za-z0-9_.:-]+", "_", x)
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(1 - eps, pmax(eps, p))
}

sigmoid_cline <- function(x, center, width, pMin, pMax) {
  pMin + (pMax - pMin) * (1 / (1 + exp(-4 * (x - center) / width)))
}

failure_row <- function(path_tsv, stage, reason) {
  tibble(
    locus_id = path_file(path_tsv) %>% path_ext_remove(),
    file = as.character(path_tsv),
    stage = stage,
    reason = reason
  )
}

read_one_locus <- function(path_tsv) {
  df <- tryCatch(
    read.table(path_tsv, header = TRUE, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)

  req <- c("site_id", "altitude", "k_minor", "n_chrom")
  if (!all(req %in% names(df))) return(NULL)

  df <- df %>%
    mutate(
      altitude = as.numeric(altitude),
      k_minor  = as.numeric(k_minor),
      n_chrom  = as.numeric(n_chrom)
    ) %>%
    filter(
      is.finite(altitude),
      is.finite(k_minor),
      is.finite(n_chrom),
      n_chrom > 0,
      k_minor >= 0,
      k_minor <= n_chrom
    ) %>%
    arrange(altitude)

  if (nrow(df) < 3L) return(NULL)

  df %>% mutate(p = clip_prob(k_minor / n_chrom))
}

trend_from_observed <- function(x, p, threshold = rho_threshold) {
  rho <- suppressWarnings(cor(x, p, method = "spearman"))
  trend <- dplyr::case_when(
    is.finite(rho) & rho >  threshold ~ "Increasing",
    is.finite(rho) & rho < -threshold ~ "Decreasing",
    TRUE                              ~ "Flat"
  )
  list(trend = trend, rho = rho)
}

make_bounds <- function(x) {
  alt_min <- min(x, na.rm = TRUE)
  alt_max <- max(x, na.rm = TRUE)
  alt_rng <- alt_max - alt_min
  list(
    center_low = alt_min - center_buffer_prop * alt_rng,
    center_high = alt_max + center_buffer_prop * alt_rng,
    width_max = width_max_factor * alt_rng
  )
}

apply_hzar_bounds <- function(model, bounds) {
  if (!use_parameter_bounds) return(model)
  if (is.null(bounds)) return(model)

  model2 <- try(
    hzar.model.addCenterRange(model, low = bounds$center_low, high = bounds$center_high),
    silent = TRUE
  )
  if (!inherits(model2, "try-error")) model <- model2

  model2 <- try(
    hzar.model.addMaxWidth(model, maxValue = bounds$width_max),
    silent = TRUE
  )
  if (!inherits(model2, "try-error")) model <- model2

  model
}

fit_hzar_basic <- function(x, p, n, burnin = main_burnin, mcmc = main_mcmc, bounds = NULL) {
  obs <- hzar.doMolecularData1DPops(distance = x, pObs = p, n)
  model <- hzar.makeCline1DFreq(obs, scaling = "free", tails = "none")
  model <- apply_hzar_bounds(model, bounds)
  req <- hzar.first.fitRequest.old.ML(model, obs)

  if (!is.null(req$mcmcParam)) {
    req$mcmcParam$burnin <- as.integer(burnin)
    req$mcmcParam$chainLength <- as.integer(mcmc)
  }

  fit <- try(hzar.doFit(req), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  ml <- try(hzar.get.ML.cline(fit), silent = TRUE)
  if (inherits(ml, "try-error")) return(NULL)

  list(fit = fit, ml = ml)
}

extract_params <- function(ml) {
  p_all <- try(ml$param.all, silent = TRUE)
  if (inherits(p_all, "try-error") || is.null(p_all)) return(NULL)

  out <- list(
    center = suppressWarnings(as.numeric(p_all$center)),
    width  = suppressWarnings(as.numeric(p_all$width)),
    pMin   = suppressWarnings(as.numeric(p_all$pMin)),
    pMax   = suppressWarnings(as.numeric(p_all$pMax))
  )

  if (any(!is.finite(unlist(out)))) return(NULL)
  out
}

empty_boot <- function(n_boot) {
  c(
    center_lo = NA_real_, center_hi = NA_real_,
    width_lo = NA_real_, width_hi = NA_real_,
    delta_p_abs_lo = NA_real_, delta_p_abs_hi = NA_real_,
    n_boot_total = as.numeric(n_boot), n_boot_ok = 0,
    boot_success_prop = NA_real_,
    center_boot_sd = NA_real_, width_boot_sd = NA_real_, delta_p_abs_boot_sd = NA_real_,
    center_ci_hits_bound = NA_real_, width_ci_hits_bound = NA_real_
  )
}

summarise_boot_draws <- function(draws, n_boot, bounds = NULL) {
  empty <- empty_boot(n_boot)
  ok <- complete.cases(draws[, c("center", "width", "delta_p_abs"), drop = FALSE])
  n_ok <- sum(ok)
  success_prop <- n_ok / n_boot

  if (n_ok == 0L) {
    empty["boot_success_prop"] <- success_prop
    return(empty)
  }

  c_q <- quantile(draws[ok, "center"], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  w_q <- quantile(draws[ok, "width"],  c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  d_q <- quantile(draws[ok, "delta_p_abs"], c(0.025, 0.975), na.rm = TRUE, names = FALSE)

  near_bound <- function(vals, low = NA_real_, high = NA_real_, tol = 1e-6) {
    out <- FALSE
    if (is.finite(low)) out <- out | any(abs(vals - low) <= tol, na.rm = TRUE)
    if (is.finite(high)) out <- out | any(abs(vals - high) <= tol, na.rm = TRUE)
    out
  }

  center_bound_hit <- if (!is.null(bounds)) near_bound(draws[ok, "center"], bounds$center_low, bounds$center_high) else NA
  width_bound_hit <- if (!is.null(bounds)) near_bound(draws[ok, "width"], NA_real_, bounds$width_max) else NA

  c(
    center_lo = c_q[[1]],
    center_hi = c_q[[2]],
    width_lo = w_q[[1]],
    width_hi = w_q[[2]],
    delta_p_abs_lo = d_q[[1]],
    delta_p_abs_hi = d_q[[2]],
    n_boot_total = as.numeric(n_boot),
    n_boot_ok = as.numeric(n_ok),
    boot_success_prop = success_prop,
    center_boot_sd = ifelse(n_ok > 1L, sd(draws[ok, "center"], na.rm = TRUE), NA_real_),
    width_boot_sd = ifelse(n_ok > 1L, sd(draws[ok, "width"], na.rm = TRUE), NA_real_),
    delta_p_abs_boot_sd = ifelse(n_ok > 1L, sd(draws[ok, "delta_p_abs"], na.rm = TRUE), NA_real_),
    center_ci_hits_bound = as.numeric(center_bound_hit),
    width_ci_hits_bound = as.numeric(width_bound_hit)
  )
}

boot_hzar_params_site <- function(x, p_fit, n, n_boot = 100L, bounds = NULL) {
  if (length(x) < 3L || n_boot <= 0L) return(empty_boot(n_boot))

  draws <- matrix(NA_real_, nrow = n_boot, ncol = 3L)
  colnames(draws) <- c("center", "width", "delta_p_abs")

  for (bi in seq_len(n_boot)) {
    idx <- sample.int(length(x), replace = TRUE)
    res <- fit_hzar_basic(x[idx], p_fit[idx], n[idx], burnin = boot_burnin, mcmc = boot_mcmc, bounds = bounds)
    if (is.null(res)) next

    pa <- extract_params(res$ml)
    if (is.null(pa)) next

    draws[bi, "center"] <- as.numeric(pa$center)
    draws[bi, "width"] <- as.numeric(pa$width)
    draws[bi, "delta_p_abs"] <- abs(as.numeric(pa$pMax - pa$pMin))
  }

  summarise_boot_draws(draws, n_boot, bounds)
}

boot_hzar_params_parametric <- function(x, p_fit, n, params, n_boot = 100L, bounds = NULL) {
  if (length(x) < 3L || n_boot <= 0L) return(empty_boot(n_boot))

  p_hat <- sigmoid_cline(x, params$center, params$width, params$pMin, params$pMax)
  p_hat <- pmin(1 - 1e-12, pmax(1e-12, p_hat))

  draws <- matrix(NA_real_, nrow = n_boot, ncol = 3L)
  colnames(draws) <- c("center", "width", "delta_p_abs")

  for (bi in seq_len(n_boot)) {
    k_boot <- stats::rbinom(length(n), size = as.integer(round(n)), prob = p_hat)
    p_boot <- clip_prob(k_boot / n)

    res <- fit_hzar_basic(x, p_boot, n, burnin = boot_burnin, mcmc = boot_mcmc, bounds = bounds)
    if (is.null(res)) next

    pa <- extract_params(res$ml)
    if (is.null(pa)) next

    draws[bi, "center"] <- as.numeric(pa$center)
    draws[bi, "width"] <- as.numeric(pa$width)
    draws[bi, "delta_p_abs"] <- abs(as.numeric(pa$pMax - pa$pMin))
  }

  summarise_boot_draws(draws, n_boot, bounds)
}

fit_one_file <- function(path_tsv) {
  df <- read_one_locus(path_tsv)
  if (is.null(df)) {
    return(list(result = NULL, failure = failure_row(path_tsv, "read_locus", "missing columns, invalid values, or fewer than three usable sites")))
  }

  x <- df$altitude
  p_orig <- df$p
  n <- df$n_chrom
  k_orig <- df$k_minor

  trend_info <- trend_from_observed(x, p_orig)
  fit_flipped <- orient_to_increasing && identical(trend_info$trend, "Decreasing")

  k_fit <- if (fit_flipped) n - k_orig else k_orig
  p_fit <- clip_prob(k_fit / n)
  bounds <- make_bounds(x)

  res_basic <- fit_hzar_basic(x, p_fit, n, bounds = bounds)
  if (is.null(res_basic)) {
    return(list(result = NULL, failure = failure_row(path_tsv, "main_fit", "hzar main fit failed")))
  }

  params <- extract_params(res_basic$ml)
  if (is.null(params)) {
    return(list(result = NULL, failure = failure_row(path_tsv, "extract_params", "could not extract finite hzar parameters")))
  }

  p_hat_fit <- sigmoid_cline(x, params$center, params$width, params$pMin, params$pMax)
  eps <- 1e-12
  p_hat_fit <- pmin(1 - eps, pmax(eps, p_hat_fit))
  ll <- sum(k_fit * log(p_hat_fit) + (n - k_fit) * log(1 - p_hat_fit))
  aic <- 2 * 4 - 2 * ll

  ci <- if (do_bootstrap && bootstrap_type == "parametric") {
    boot_hzar_params_parametric(x, p_fit, n, params, n_boot = n_boot, bounds = bounds)
  } else if (do_bootstrap && bootstrap_type == "site") {
    boot_hzar_params_site(x, p_fit, n, n_boot = n_boot, bounds = bounds)
  } else {
    empty_boot(0L)
  }

  p_low_orig <- if (fit_flipped) 1 - as.numeric(params$pMin) else as.numeric(params$pMin)
  p_high_orig <- if (fit_flipped) 1 - as.numeric(params$pMax) else as.numeric(params$pMax)
  delta_p_signed <- p_high_orig - p_low_orig

  out <- tibble(
    locus_id = path_file(path_tsv) %>% path_ext_remove(),
    model = "free_noTails",
    bootstrap_type = ifelse(do_bootstrap, bootstrap_type, "none"),
    fit_orientation = ifelse(fit_flipped, "flipped_to_increasing", "original"),
    trend_spearman = trend_info$trend,
    rho_spearman = as.numeric(trend_info$rho),
    n_sites = nrow(df),
    alt_min = min(x, na.rm = TRUE),
    alt_max = max(x, na.rm = TRUE),
    alt_rng = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
    center_bound_low = bounds$center_low,
    center_bound_high = bounds$center_high,
    width_bound_high = bounds$width_max,
    dp_obs = diff(range(k_orig / n, na.rm = TRUE)),
    center = as.numeric(params$center),
    center_lo = as.numeric(ci["center_lo"]),
    center_hi = as.numeric(ci["center_hi"]),
    width = as.numeric(params$width),
    width_lo = as.numeric(ci["width_lo"]),
    width_hi = as.numeric(ci["width_hi"]),
    pMin = as.numeric(params$pMin),
    pMax = as.numeric(params$pMax),
    p_low_orig = p_low_orig,
    p_high_orig = p_high_orig,
    delta_p = as.numeric(params$pMax - params$pMin),
    delta_p_signed = delta_p_signed,
    delta_p_abs = abs(delta_p_signed),
    delta_p_abs_lo = as.numeric(ci["delta_p_abs_lo"]),
    delta_p_abs_hi = as.numeric(ci["delta_p_abs_hi"]),
    logLik = ll,
    AIC = aic,
    n_boot_total = as.numeric(ci["n_boot_total"]),
    n_boot_ok = as.numeric(ci["n_boot_ok"]),
    boot_success_prop = as.numeric(ci["boot_success_prop"]),
    center_boot_sd = as.numeric(ci["center_boot_sd"]),
    width_boot_sd = as.numeric(ci["width_boot_sd"]),
    delta_p_abs_boot_sd = as.numeric(ci["delta_p_abs_boot_sd"]),
    center_ci_hits_bound = as.logical(as.numeric(ci["center_ci_hits_bound"])),
    width_ci_hits_bound = as.logical(as.numeric(ci["width_ci_hits_bound"]))
  )

  list(result = out, failure = NULL)
}

## ====== step 1: generate per-locus hzar input files =================================
if (force_rebuild_inputs && dir_exists(out_hzar_dir)) {
  message("force rebuilding per-locus hzar input files: deleting ", out_hzar_dir)
  dir_delete(out_hzar_dir)
}

if (!dir_exists(out_hzar_dir) || length(dir_ls(out_hzar_dir)) == 0L) {
  message("Generating per-locus HZAR input files...")

  geno_raw <- read_tsv(
    geno_path,
    col_types = cols(
      Chr = col_character(),
      Pos = col_character(),
      Maj = col_character(),
      Min = col_character(),
      .default = col_character()
    ),
    trim_ws = TRUE,
    na = c("NA", "Na", "na", ".", "")
  )

  info <- read_tsv(info_path, col_types = cols(.default = col_character()))

  stopifnot(all(c("Chr", "Pos", "Maj", "Min") %in% names(geno_raw)))
  stopifnot(all(c("sample_id", "altitude") %in% names(info)))

  info <- info %>% mutate(altitude = as.numeric(altitude))
  if (!"site_id" %in% names(info)) info <- info %>% mutate(site_id = paste0("ALT_", altitude))

  sample_cols <- setdiff(names(geno_raw), c("Chr", "Pos", "Maj", "Min"))
  stopifnot(length(sample_cols) > 0L)

  geno_long <- geno_raw %>%
    pivot_longer(all_of(sample_cols), names_to = "sample_id", values_to = "GT") %>%
    mutate(
      Chr = as.character(Chr),
      Pos = as.character(Pos),
      locus_id = paste(Chr, Pos, sep = "_"),
      Maj = toupper(Maj),
      Min = toupper(Min)
    ) %>%
    inner_join(info[, c("sample_id", "altitude", "site_id")], by = "sample_id")

  allele_check <- geno_long %>%
    group_by(locus_id) %>%
    summarise(nMaj = n_distinct(Maj), nMin = n_distinct(Min), .groups = "drop") %>%
    filter(nMaj != 1L | nMin != 1L)

  if (nrow(allele_check) > 0L) {
    warning(
      "loci with inconsistent Maj/Min will be skipped: ",
      paste(head(allele_check$locus_id, 10), collapse = ", "),
      if (nrow(allele_check) > 10L) " ..." else ""
    )
  }

  dir_create(out_hzar_dir)
  input_log <- vector("list", length(unique(geno_long$locus_id)))
  n_written <- 0L
  idx_log <- 1L

  for (loc in unique(geno_long$locus_id)) {
    sub <- geno_long %>% filter(locus_id == loc)
    maj_allele <- sub %>% distinct(Maj) %>% pull()
    min_allele <- sub %>% distinct(Min) %>% pull()

    if (length(maj_allele) != 1L || length(min_allele) != 1L) {
      input_log[[idx_log]] <- tibble(locus_id = loc, status = "skipped", reason = "inconsistent Maj/Min", n_sites = NA_integer_, n_chrom_total = NA_integer_, global_minor_p = NA_real_, recoded = NA)
      idx_log <- idx_log + 1L
      next
    }

    g_counts <- convert_to_minor_counts(sub$GT, maj = maj_allele[[1]], min = min_allele[[1]])
    n_called_genotypes <- sum(!is.na(g_counts))

    if (n_called_genotypes == 0L) {
      input_log[[idx_log]] <- tibble(locus_id = loc, status = "skipped", reason = "no called genotypes", n_sites = NA_integer_, n_chrom_total = 0L, global_minor_p = NA_real_, recoded = NA)
      idx_log <- idx_log + 1L
      next
    }

    alt_k_global <- sum(g_counts, na.rm = TRUE)
    n_chr_global <- 2L * n_called_genotypes
    alt_p_global <- alt_k_global / n_chr_global
    recoded <- FALSE

    if (recode_to_global_minor && is.finite(alt_p_global) && alt_p_global > 0.5) {
      g_counts <- ifelse(is.na(g_counts), NA_integer_, 2L - g_counts)
      recoded <- TRUE
      alt_p_global <- 1 - alt_p_global
    }

    sub$k_minor <- g_counts
    sub$n_chrom <- ifelse(is.na(g_counts), 0L, 2L)

    per_site <- sub %>%
      filter(!is.na(altitude)) %>%
      group_by(site_id, altitude) %>%
      summarise(
        k_minor = sum(k_minor, na.rm = TRUE),
        n_chrom = sum(n_chrom, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(n_chrom > 0) %>%
      arrange(altitude)

    n_sites <- nrow(per_site)
    if (n_sites < min_sites_per_locus) {
      input_log[[idx_log]] <- tibble(locus_id = loc, status = "skipped", reason = "fewer than min_sites_per_locus usable sites", n_sites = n_sites, n_chrom_total = sum(per_site$n_chrom), global_minor_p = alt_p_global, recoded = recoded)
      idx_log <- idx_log + 1L
      next
    }

    out_file <- fs::path(out_hzar_dir, paste0(sanitize_locus_id(loc), ".tsv"))
    write_tsv_plain(per_site, out_file)
    n_written <- n_written + 1L

    input_log[[idx_log]] <- tibble(locus_id = loc, status = "written", reason = NA_character_, n_sites = n_sites, n_chrom_total = sum(per_site$n_chrom), global_minor_p = alt_p_global, recoded = recoded)
    idx_log <- idx_log + 1L
  }

  input_log_tbl <- bind_rows(input_log[!vapply(input_log, is.null, logical(1))])
  write_tsv_plain(input_log_tbl, out_input_log_tsv)
  message("Built ", n_written, " per-locus TSVs in ", out_hzar_dir)
  message("Wrote: ", out_input_log_tsv)
} else {
  message("Using existing per-locus HZAR input files in ", out_hzar_dir)
  message("Use --force_rebuild_inputs to regenerate them from genotype and metadata files.")
}

## ====== step 2: choose loci to fit ===================================================
all_files <- sort(dir_ls(out_hzar_dir, glob = "*.tsv"))
stopifnot(length(all_files) > 0L)

if (!is.na(max_loci)) {
  all_files <- head(all_files, max_loci)
  message("Testing mode: limiting run to first ", length(all_files), " loci because --max_loci was provided.")
}

tsv_files <- all_files
current_out_tsv <- out_param_tsv
current_fail_tsv <- "hzar_fit_failures.tsv"
message("Processing all ", length(tsv_files), " loci and writing one output table: ", current_out_tsv)

## ====== step 3: execute hzar fitting ================================================
fit_records <- mclapply(tsv_files, fit_one_file, mc.cores = n_cores, mc.preschedule = TRUE, mc.set.seed = TRUE)

result_list <- lapply(fit_records, function(x) x$result)
failure_list <- lapply(fit_records, function(x) x$failure)

result_list <- result_list[!vapply(result_list, is.null, logical(1))]
failure_list <- failure_list[!vapply(failure_list, is.null, logical(1))]

if (length(result_list) > 0L) {
  res <- bind_rows(result_list)
  write_tsv_plain(res %>% arrange(locus_id), current_out_tsv)
  message("Cline fitting complete. Wrote file: ", current_out_tsv)
} else {
  message("No models successfully fitted in this run.")
}

if (length(failure_list) > 0L) {
  failures <- bind_rows(failure_list)
  write_tsv_plain(failures %>% arrange(locus_id), current_fail_tsv)
  message("Wrote fit failure log: ", current_fail_tsv)
}
