# Authors: Ufuk Topalan, Ismail Kudret Saglam
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

setwd(".")

# ---- Load likelihood table ----
df <- read.table("admix_runs_LH.txt", header = TRUE)

# Summaries per K
sumK <- df %>%
  group_by(K) %>%
  summarise(
    mean_LH = mean(LH),
    sd_LH   = sd(LH),
    var_LH  = var(LH),
    n = n()
  )

print(sumK)

sumK <- sumK %>% arrange(K)
Kmax <- max(sumK$K)

# compute Evanno following the symmetric second-difference form:
# L''(K) = L(K+1) - 2*L(K) + L(K-1)
# ΔK is defined for K = 2, ..., Kmax-1

sumK2 <- sumK %>%
  mutate(
    mean_LH_lag  = lag(mean_LH),   # L(K-1)
    mean_LH_lead = lead(mean_LH)   # L(K+1)
  ) %>%
  mutate(
    L_doubleprime = mean_LH_lead - 2*mean_LH + mean_LH_lag
  ) %>%
  # Evanno: deltaK = |L''(K)| / sd(L(K))
  # but only meaningful for K = 2 .. Kmax-1
  mutate(
    deltaK_raw = ifelse(!is.na(L_doubleprime) & sd_LH > 0,
                        abs(L_doubleprime) / sd_LH,
                        NA_real_)
  ) %>%
  select(K, mean_LH, sd_LH, var_LH, n, L_doubleprime, deltaK = deltaK_raw)

print(sumK2)

# LnP(K) plot (with sd bars)
p1 <- ggplot(sumK, aes(x = K, y = mean_LH)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_LH - sd_LH, ymax = mean_LH + sd_LH), width = 0.15) +
  theme_minimal(base_size = 14) +
  labs(title = "Mean LnP(K) with SD", y = "Mean LnP(K)")

# ΔK plot: only plot K where deltaK is finite
# Choose your preferred color:
mycol <- "#009E73"   # green (colorblind-friendly, no yellow)

# Filter only finite deltaK values for labeling
df_labels <- df_plot %>%
  filter(is.finite(deltaK)) %>%
  mutate(label = sprintf("%.2f", deltaK))   # 2 decimals; switch to "%.3f" if you want 3 decimals

p2 <- ggplot(df_plot, aes(x = K, y = deltaK)) +
  geom_line(color = mycol, size = 1, na.rm = TRUE) +
  geom_point(color = mycol, size = 3, na.rm = TRUE) +
  
  # Add text labels for deltaK values
  geom_text(
    data = df_labels,
    aes(label = label),
    vjust = -0.6,       # moves label slightly above the point
    size = 4,           # text size (adjust if needed)
    color = "black"
  ) +
  
  scale_x_continuous(
    breaks = Kmin:Kmax,
    limits = c(Kmin, Kmax)
  ) +
  labs(
    title = "Evanno ΔK",
    x = "K",
    y = "ΔK"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )


print(p1)
print(p2)

candidates <- sumK2 %>%
  filter(!is.na(deltaK)) %>%
  arrange(desc(deltaK))

print(candidates)

bestK_by_evanno <- candidates$K[which.max(candidates$deltaK)]
cat("Best K by Evanno (where defined):", bestK_by_evanno, "\n")
