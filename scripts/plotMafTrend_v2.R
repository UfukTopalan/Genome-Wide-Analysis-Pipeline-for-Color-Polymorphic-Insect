############################# Allele Frequency Trend #############################

# Input (minor allele frequencies across altitudes)
freq_path <- "isophya71.alt.assoc.loci.mafs.tsv"

# Output
out_meta_pdf <- "isophya71.alt.assoc.maf.trend_v3.pdf"

## ====== LIBRARIES =============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

## --------------------------------------------------------------------------
## PART 1. META-CLINE (two-slope regression)
## --------------------------------------------------------------------------

freq <- read.table(freq_path,
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

stopifnot(all(c("chromo", "position") %in% names(freq)))

# Extract altitude columns only
alt_map <- function(df) {

  cols <- names(df)

  tibble(
    colname = cols,
    altitude = suppressWarnings(as.numeric(cols))
  ) %>%
    filter(!is.na(altitude))
}

freq_cols <- alt_map(freq)

freq_long <- freq %>%
  pivot_longer(
    all_of(freq_cols$colname),
    names_to = "alt_col",
    values_to = "p"
  ) %>%
  left_join(freq_cols,
            by = c("alt_col" = "colname")) %>%
  transmute(
    chromo,
    position,
    altitude = altitude,
    p = as.numeric(p)
  ) %>%
  filter(is.finite(p)) %>%
  mutate(
    p = pmin(1, pmax(0, p))
  )

# --------------------------------------------------------------------------
# Classify loci by overall allele-frequency trend
# --------------------------------------------------------------------------

trend_by_locus <- freq_long %>%
  group_by(chromo, position) %>%
  summarise(
    fit = list(lm(p ~ altitude)),
    .groups = "drop"
  ) %>%
  mutate(
    slope = purrr::map_dbl(fit, ~coef(.x)[2]),
    trend = case_when(
      slope > 1e-6  ~ "Increasing",
      slope < -1e-6 ~ "Decreasing",
      TRUE ~"Flat"
    )
  ) %>%
  select(-fit)

# One observation per locus per altitude
locus_means <- freq_long %>%
  inner_join(
    trend_by_locus,
    by = c("chromo", "position")
  ) %>%
  rename(
    p_locus = p
  )

# --------------------------------------------------------------------------
# Test whether increasing and decreasing loci have different slopes
# --------------------------------------------------------------------------

fit_meta <- lm(
  p_locus ~ altitude * trend,
  data = locus_means
)

print(summary(fit_meta))

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

summary_by_alt <- locus_means %>%
  group_by(trend, altitude) %>%
  summarise(
    mean_p = mean(p_locus),
    se = sd(p_locus) / sqrt(n()),
    .groups = "drop"
  )

p_two_slope <- ggplot(
  summary_by_alt,
  aes(
    x = altitude,
    y = mean_p,
    color = trend,
    fill = trend
  )
) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(
    aes(
      ymin = mean_p - 1.96 * se,
      ymax = mean_p + 1.96 * se
    ),
    alpha = 0.20,
    color = NA
  ) +
  scale_color_manual(
    values = c(
      Increasing = "#00441B",
      Decreasing = "#66C2A5"
    )
  ) +
  scale_fill_manual(
    values = c(
      Increasing = "#00441B",
      Decreasing = "#66C2A5"
    )
  ) +
  labs(
    x = "Altitude (m)",
    y = "Mean frequency of globally defined minor allele",
    color = "Locus group",
    fill = "Locus group"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  out_meta_pdf,
  p_two_slope,
  width = 140,
  height = 90,
  units = "mm"
)

message("Wrote: ", out_meta_pdf)
