# Authors: Ufuk Topalan, Ismail Kudret Saglam
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------
## admixture_plot_isophya.R


setwd(".")

library(tidyverse)
library(cowplot)

## 1. metadata -------------------------------------------------------------

info <- read.table("isophya71.info", header = FALSE, stringsAsFactors = FALSE)
colnames(info) <- c("pop", "id")

## fixed altitude order
alt_order <- c(450, 850, 900, 1000, 1100, 1200, 1300, 1900, 2000, 2100, 2300)
info$pop <- factor(info$pop, levels = alt_order)

## base palette for "other" clusters (beyond the low/high pair)
other_palette <- c(
  "#66C2A4",  # teal
  "#8DA0CB",  # bluish
  "#B2DF8A"   # light green
)

## 2. helper: read Q, enforce order, clamp to [0,1], define colors --------

make_q_df <- function(qfile, info) {
  qmat <- as.matrix(read.table(qfile))
  K <- ncol(qmat)
  colnames(qmat) <- paste0("C", seq_len(K))
  
  ## clamp any stray values to [0, 1]
  qmat[qmat < 0] <- 0
  qmat[qmat > 1] <- 1
  
  ## replace NA with 0 so geom_col() never drops rows
  qmat[is.na(qmat)] <- 0
  
  ## renormalize rows to sum 1 (protect against rounding / NA fix)
  rs <- rowSums(qmat)
  ok <- which(rs > 0)
  qmat[ok, ] <- sweep(qmat[ok, , drop = FALSE], 1, rs[ok], "/")
  
  df_wide <- cbind(info, qmat) %>%
    arrange(pop, id) %>%
    mutate(ind_order = row_number())
  
  bounds <- df_wide %>%
    group_by(pop) %>%
    summarise(
      min_ind = min(ind_order),
      max_ind = max(ind_order),
      mid_ind = (min_ind + max_ind) / 2,
      .groups = "drop"
    )
  
  separators <- bounds$max_ind[-nrow(bounds)] + 0.5
  
  df_long <- df_wide %>%
    pivot_longer(starts_with("C"),
                 names_to = "cluster",
                 values_to = "q")
  
  df_long$cluster <- factor(df_long$cluster, levels = paste0("C", seq_len(K)))
  
  low_pops  <- c(450, 850, 900, 1000)
  high_pops <- setdiff(alt_order, low_pops)
  
  low_df  <- df_wide %>% filter(pop %in% low_pops)
  high_df <- df_wide %>% filter(pop %in% high_pops)
  
  cluster_names <- paste0("C", seq_len(K))
  means_low  <- colMeans(low_df[, cluster_names, drop = FALSE])
  means_high <- colMeans(high_df[, cluster_names, drop = FALSE])
  
  cluster_low  <- names(which.max(means_low))
  cluster_high <- names(which.max(means_high))
  
  color_map <- setNames(rep(NA_character_, K), cluster_names)
  color_map[cluster_low]  <- "#000000"
  color_map[cluster_high] <- "#009E73"
  
  other_clusters <- setdiff(cluster_names, c(cluster_low, cluster_high))
  if (length(other_clusters) > 0) {
    color_map[other_clusters] <- other_palette[seq_along(other_clusters)]
  }
  
  list(
    df_long    = df_long,
    df_wide    = df_wide,
    bounds     = bounds,
    separators = separators,
    K          = K,
    color_map  = color_map
  )
}


## 3. Figure 2B: K = 2 main panel -----------------------------------------

plot_admix_k2_main <- function(qfile, info) {
  x <- make_q_df(qfile, info)
  d <- x$df_long
  bounds <- x$bounds
  separators <- x$separators
  color_map <- x$color_map
  
  p <- ggplot(d, aes(x = ind_order, y = q, fill = cluster)) +
    geom_col(width = 1) +
    geom_vline(
      xintercept = separators,
      linetype   = "dashed",
      linewidth  = 0.3,
      colour     = "white"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = bounds$mid_ind,
      labels = as.character(bounds$pop)
    ) +
    scale_fill_manual(values = color_map, guide = "none") +  # no legend
    labs(
      x = "Altitude (m)",
      y = "Admixture proportions"
    ) +
    theme_bw() +
    theme(
      panel.grid    = element_blank(),
      axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 7),
      axis.text.y   = element_text(size = 8),
      axis.title    = element_text(size = 9),
      legend.position = "none",
      plot.margin   = margin(t = 5, r = 5, b = 10, l = 5)
    )
  
  p
}

p_k2 <- plot_admix_k2_main("isophya71.admix2.run1.qopt", info)
ggsave("isophya71_admixture_K2.pdf", p_k2, width = 6, height = 4)


## 4. Multi-K figure: K = 2–5 stacked ------------------------------------

plot_admix_panel <- function(qfile, info, title, show_ylab = FALSE) {
  x <- make_q_df(qfile, info)
  d <- x$df_long
  bounds <- x$bounds
  separators <- x$separators
  color_map <- x$color_map
  
  p <- ggplot(d, aes(x = ind_order, y = q, fill = cluster)) +
    geom_col(width = 1) +
    geom_vline(
      xintercept = separators,
      linetype   = "dashed",
      linewidth  = 0.3,
      colour     = "white"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = bounds$mid_ind,
      labels = as.character(bounds$pop)
    ) +
    scale_fill_manual(values = color_map, guide = "none") +
    labs(
      x = "Sampling altitude (m a.s.l.)",
      y = if (show_ylab) "Admixture proportions" else NULL,
      title = title
    ) +
    theme_bw() +
    theme(
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0, size = 10, face = "bold"),
      axis.text.x  = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8),
      axis.text.y  = element_text(size = 8.5),
      axis.title.x = element_text(size = 9.5),
      axis.title.y = if (show_ylab) {
        element_text(size = 9.5)
      } else {
        element_text(colour = NA, size = 9.5)
      },
      legend.position = "none",
      plot.margin  = margin(t = 6, r = 5, b = 6, l = 5)
    )
  
  p
}

p2 <- plot_admix_panel("isophya71.admix2.run1.qopt",
                       info, title = "K = 2", show_ylab = TRUE)
p3 <- plot_admix_panel("isophya71_admix3_run1.qopt",
                       info, title = "K = 3", show_ylab = FALSE)
p4 <- plot_admix_panel("isophya71_admix4_run1.qopt",
                       info, title = "K = 4", show_ylab = FALSE)
p5 <- plot_admix_panel("isophya71_admix5_run1.qopt",
                       info, title = "K = 5", show_ylab = FALSE)

p_multi <- plot_grid(p2, p3, p4, p5,
                     ncol = 1,
                     align = "v",
                     rel_heights = c(1, 1, 1, 1))

ggsave("isophya71_admixture_K2_K5.pdf", p_multi, width = 7, height = 4.5)
