# plot_all_folded_sfs.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# normalize helper
norm <- function(x) x / sum(x)

# read + clean one folded SFS file
read_folded_sfs <- function(path) {
  sfs <- scan(path, quiet = TRUE)

  if (length(sfs) < 3) stop("SFS too short in: ", path)

  # drop the fixed bins (your original behavior): first and last
  sfs_var <- sfs[-c(1, length(sfs))]

  # trim trailing zeros that come from unfolded-length padding
  # keep everything up to the last non-zero entry
  nz <- which(sfs_var != 0)
  if (length(nz) == 0) stop("All variable bins are zero after trimming fixed bins in: ", path)
  sfs_var <- sfs_var[seq_len(max(nz))]

  tibble(
    pop = str_remove(basename(path), "\\.sfs$") |> str_remove("_folded$"),
    maf_bin = seq_along(sfs_var),
    prop = norm(sfs_var)
  )
}

# --- input: all *_folded.sfs files in a directory ---
sfs_dir <- "."
files <- list.files(sfs_dir, pattern = "_folded\\.sfs$", full.names = TRUE)

stopifnot(length(files) > 0)

df <- bind_rows(lapply(files, read_folded_sfs))

# order facets
df <- df %>%
  mutate(order_key = suppressWarnings(as.numeric(str_extract(pop, "\\d+")))) %>%
  arrange(order_key, pop) %>%
  mutate(pop = factor(pop, levels = unique(pop)))

p <- ggplot(df, aes(x = maf_bin, y = prop)) +
  geom_col(fill = "#00441B") +
  facet_wrap(~ pop, scales = "free_x") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1),
    minor_breaks = NULL
  ) +
  labs(x = "number of minor alleles", y = "Proportions") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title  = element_text(size = 11)
  )

print(p)

ggsave(
  "SFS_plots_faceted_2over3_A4.png",
  p,
  width  = 210,
  height = 198,
  units  = "mm",
  dpi    = 300
)

ggsave(
  "SFS_plots_faceted_2over3_A4.pdf",
  p,
  width  = 210,
  height = 198,
  units  = "mm",
  dpi    = 300
)
