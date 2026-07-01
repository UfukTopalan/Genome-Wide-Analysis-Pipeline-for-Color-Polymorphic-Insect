library(vegan)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# Load the data
neutral <- read.table("isophya71.neut.fst.tsv", header = TRUE)
adaptive <- read.table("isophya71.adapt.fst.tsv", header = TRUE)

# Calculate linearized Fst for both datasets
neutral$Linearized_Fst <- neutral$Fst / (1 - neutral$Fst)
adaptive$Linearized_Fst <- adaptive$Fst / (1 - adaptive$Fst)

# Create regression plot with separate regression lines for adaptive and neutral loci
neutral$Loci <- "Neutral"
adaptive$Loci <- "Adaptive"
combined <- rbind(neutral, adaptive)

p <- ggplot(combined, aes(x = Distance, y = Linearized_Fst, color = Loci)) +
  geom_point(alpha = 0.7, size = 3) +  # Make the dots bigger
  geom_smooth(data = subset(combined, Loci == "Adaptive"), method = "lm", se = FALSE, color = "#006D2C") +  # Adaptive regression line
  geom_smooth(data = subset(combined, Loci == "Neutral"), method = "lm", se = FALSE, color = "#66C2A5") +  # Neutral regression line
  scale_color_manual(values = c("Adaptive" = "#006D2C", "Neutral" = "#66C2A5")) +  # Custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),  # Adjust x-axis text size
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis label size
    axis.title.y = element_text(size = 12),  # Adjust y-axis label size
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10),   # Adjust legend text size
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(size = 0.7, color = "black"),  # Make axes visible
    axis.ticks = element_line(size = 0.7, color = "black"),
    plot.margin = margin(5, 5, 5, 5)  # Adjust margins for visibility
  ) +
  labs(
    x = "Geographic Distance (m)",  # Adding unit for clarity
    y = "Fst/(1-Fst)",
    color = "Loci Type"
  )+
  theme(
    aspect.ratio = 0.6  # Adjust aspect ratio
  )

# Save the plot
ggsave("isophya71.neut.adapt.ibd.pdf", plot = p, width = 12, height = 8, units = "cm", dpi = 300)
