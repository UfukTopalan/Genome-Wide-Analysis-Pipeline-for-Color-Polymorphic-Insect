# Load necessary libraries
library(ggplot2)

# Load the data
data <- read.table("isophya71.fst.dist.tsv2", header = TRUE)
data$Linearized_Fst <- data$Fst / (1 - data$Fst)


# Create the plot
plot <- ggplot(data, aes(x = Distance, y = Linearized_Fst)) +
  geom_point(alpha = 0.6, size = 3, color = "#006D2C") +  # Larger points for better visibility
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Regression line in colorblind-friendly blue
  theme_minimal(base_size = 10) +  # Minimal theme with base font size
  theme(
    axis.line = element_line(color = "black"),  # Add visible axes
    panel.grid = element_blank(),  # Remove gridlines inside the plot
    axis.text = element_text(size = 10),  # Adjust axis text size
    axis.title = element_text(size = 12),  # Adjust axis title size
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust margins
  ) +
  labs(
    x = "Distance (m)",  # X-axis label
    y = "Fst/(1-Fst)"  # Y-axis label
  )

# Display the plot
#print(plot)

# Save the plot
ggsave("isophya71.ibd.plot.pdf", plot, width = 12, height = 8, units = "cm")

