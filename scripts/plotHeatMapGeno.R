# Authors: Ufuk Topalan, Ismail Kudret Saglam
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------
# Load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the data
data <- read.table(".", header = TRUE, stringsAsFactors = FALSE)

# Transform the data into long format for plotting
data_long <- data %>%
  pivot_longer(cols = starts_with("SNP"), names_to = "SNP", values_to = "Genotype")

# Ensure individuals maintain the original order
data_long$Indv <- factor(data_long$Indv, levels = unique(data$Indv))

# Add a grouping variable for altitude segments
data_long <- data_long %>%
  mutate(Altitude = factor(Altitude, levels = unique(Altitude)))

# Calculate x-axis positions for dashed lines and labels
altitude_positions <- data_long %>%
  group_by(Altitude) %>%
  summarize(start = min(as.numeric(Indv)),
            end = max(as.numeric(Indv))) %>%
  mutate(midpoint = (start + end) / 2)

# Create a custom color palette
custom_palette <- c("0" = "#009E73", "1" = "#00441B", "2" = "#000000")

# Plot the data
plot <- ggplot(data_long, aes(x = as.numeric(Indv), y = SNP, fill = factor(Genotype))) +
  geom_tile(color = "darkgrey") +
  scale_fill_manual(values = custom_palette, name = "Genotype",
                    labels = c("Homozygous Major", "Heterozygous", "Homozygous Minor")) +
  labs(x = "Altitude", y = "Loci", title = "Genotypic states at major effect loci") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10),           # Smaller legend title
        legend.text = element_text(size = 8),            # Smaller legend text
        legend.key.size = unit(0.5, "cm"),               # Smaller legend key size
        legend.spacing.y = unit(0.1, "cm")) +            # Reduce spacing between legend items
  scale_x_continuous(breaks = altitude_positions$midpoint, 
                     labels = altitude_positions$Altitude) +
  # Add dashed vertical lines to separate altitude groups
  geom_vline(data = altitude_positions, aes(xintercept = start - 0.5), 
             linetype = "dashed", color = "white", size = 0.5) +
  # Adjust the aspect ratio to shorten the y-axis
  coord_fixed(ratio = 2.5)

# Save the plot
ggsave("genotypic_heatmap.pdf", plot, width = 12, height = 6)
