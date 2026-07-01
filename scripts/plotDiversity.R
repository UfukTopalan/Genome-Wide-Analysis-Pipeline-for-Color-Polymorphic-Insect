# Authors: Ufuk Topalan, Ismail Kudret Saglam
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------
library(ggplot2)
library(ggsignif)
library(psych)
library(pastecs)
library(viridis)
library(reshape)
library(scales)
library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(purrr)

# Load the data
Data <- read.table("isophya71.diversity.tsv", header = TRUE, sep = '\t')
Data$Type <- as.factor(Data$Type)
Data$Altitude <- as.factor(Data$Altitude)
Data[,c(5,6)] <- Data[,c(5,6)] / Data$nSites
Data <- Data[Data$nSites > 99,]

# Perform t-tests for each altitude and calculate statistics

num_comparisons <- length(unique(Data$Altitude))  # Number of unique altitude groups for correction

stats_results <- Data %>%
  group_by(Altitude) %>%
  summarise(
    t_test_result = list(t.test(tP[Type == "adaptive"], tP[Type == "neutral"])),  # Store t-test results
    mean_adaptive = mean(tP[Type == "adaptive"], na.rm = TRUE),
    median_adaptive = median(tP[Type == "adaptive"], na.rm = TRUE),
    mean_neutral = mean(tP[Type == "neutral"], na.rm = TRUE),
    median_neutral = median(tP[Type == "neutral"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Extract t-statistic, p-value, and degrees of freedom from t-test results
    t_statistic = map_dbl(t_test_result, ~ .x$statistic),
    p_value = map_dbl(t_test_result, ~ .x$p.value),
    degrees_of_freedom = map_dbl(t_test_result, ~ .x$parameter),
    bonferroni_p = p.adjust(p_value, method = "bonferroni"),  # Adjust p-values
    significance = case_when(
      bonferroni_p < 0.001 ~ "***",
      bonferroni_p < 0.01  ~ "**",
      bonferroni_p < 0.05  ~ "*",
      TRUE                 ~ NA_character_  # Replace "ns" with NA
    ),
    y_position = 0.16  # Position above boxplots
  ) %>%
  select(-t_test_result)  # Drop the t-test result column if no longer needed

write_tsv(stats_results, "isophya71.diversity.stats.tsv", quote = "none")

# Plot
p = ggplot(Data, aes(x = Altitude, y = tP, fill = Type)) +
  geom_boxplot(
    width = 0.4, notch = TRUE,
    position = position_dodge(width = 0.6)  # Increase separation
  ) +
  scale_fill_manual(values = c("adaptive" = "#006D2C", "neutral" = "#66C2A5")) +
  labs(x = "Population", y = "Theta Pi", fill = "Variation Type") +
  scale_y_continuous(
    limits = c(0, 0.18),  # Restrict y-axis range
    breaks = seq(0, 0.18, by = 0.02)  # Set step size to 0.04
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray80", size = 0.2),
    panel.grid.minor = element_line(color = "gray90", size = 0.1),
    plot.margin = margin(4, 4, 4, 4)  # Adjust margins if necessary
  ) +
  geom_text(
    data = stats_results %>% filter(!is.na(significance)),  # Exclude "ns"
    aes(x = Altitude, y = y_position, label = significance),
    inherit.aes = FALSE,
    size = 5, color = "red",
    position = position_dodge(width = 0.6)  # Align with boxplots
  ) +
  coord_flip() +  # Flip axes
  theme(
    aspect.ratio = 1.0  # Adjust aspect ratio to make the plot more compact
  )

ggsave("isophya71.diversity.pdf", plot = p, width = 18, height = 10, units = "cm", dpi = 300)

