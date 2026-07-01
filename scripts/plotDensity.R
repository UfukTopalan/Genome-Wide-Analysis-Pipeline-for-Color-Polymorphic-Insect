# Author(s): Ufuk Topalan, [Supervisor Name]
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------
library(ggplot2)

# Load the data
data <- read.csv(".")

# Create the density plot
density_plot <- ggplot(data, aes(x = Discriminant_Function_1, fill = Color)) +
  geom_density(alpha = 0.8) +
  scale_fill_manual(values = c("Dark" = "#000000", "Pale" = "#009E73")) +
  labs(x = "Discriminant Function 1", y = "Density") +
  theme_minimal(base_size = 10) + # Base font size for readability
  theme(
    axis.text = element_text(size = 9),        # Axis tick labels
    axis.title = element_text(size = 10),      # Axis titles made slightly larger for emphasis
    axis.line = element_line(size = 0.5),      # Make axis lines thicker for a traditional look
    axis.ticks = element_line(size = 0.8),     # Thicker ticks
    legend.title = element_blank(),            # Remove legend title
    legend.text = element_text(size = 8),      # Legend text
    legend.position = "top",                   # Position legend on top
    plot.margin = margin(5, 5, 5, 5),           # Reduce plot margins
    panel.grid.major = element_line(color = "gray80", size = 0.2),  # Make major gridlines lighter
    panel.grid.minor = element_line(color = "gray90", size = 0.1)   # Make minor gridlines lighter
  ) +
  xlim(-5, 5) +  # Set the x-axis limits from -5 to 5
  coord_fixed(ratio = 8) +  # Adjust aspect ratio to make the plot taller
  scale_y_continuous(expand = c(0, 0))  # Remove padding on y-axis to make it start at y = 0

# Save the plot with specific dimensions
ggsave("/home/ismail/Research/isophya/DAPC_density_plot.pdf", density_plot, width = 3.35, height = 6.70, units = "in")




















# Create the density plot
density_plot <- ggplot(data, aes(x = Discriminant_Function_1, fill = Color)) +
  geom_density(alpha = 0.8) +
  scale_fill_manual(values = c("Dark" = "#000000", "Pale" = "#009E73")) +
  labs(x = "Discriminant Function 1", y = "Density") +
  theme_minimal(base_size = 10) + # Base font size for readability
  theme(
    axis.text = element_text(size = 8),        # Axis tick labels
    axis.title = element_text(size = 10),      # Axis titles
    legend.title = element_blank(),            # Remove legend title
    legend.text = element_text(size = 8),      # Legend text
    legend.position = "top",                   # Position legend on top
    plot.margin = margin(5, 5, 5, 5)           # Reduce plot margins
  )

# Save the plot with specific dimensions
ggsave("density_plot.pdf", density_plot, width = 8.2 / 2.54, height = 6 / 2.54, units = "in")
