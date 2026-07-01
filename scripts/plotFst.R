# Authors: Ufuk Topalan, Ismail Kudret Saglam
# Written for: Genome-Wide Analysis Reveals Altitude-Associated
#              Divergence in a Color-Polymorphic Insect (Topalan & Sağlam, 2026)
# Repository: https://github.com/UfukTopalan/Genome-Wide-Analysis-Pipeline-for-Color-Polymorphic-Insect
# ------------------------------------------------------------------
library(ggplot2)
library(dplyr)

# Read data
Fst <- read.table("isophya71.fst.dist.tsv", header = TRUE, sep = "\t")

# Define altitude order (ascending: 450 to 2300)
altitude_order <- c("450", "850", "900", "1000", "1100", "1200", "1300", "1900", "2000", "2100", "2300")

# Reorder pairs to ensure Pop1 < Pop2 (lower triangle)
Fst_lower <- Fst %>%
  rowwise() %>%
  mutate(
    Pop1_new = ifelse(as.numeric(Pop1) < as.numeric(Pop2), as.character(Pop1), as.character(Pop2)),
    Pop2_new = ifelse(as.numeric(Pop1) < as.numeric(Pop2), as.character(Pop2), as.character(Pop1))
  ) %>%
  select(Pop1 = Pop1_new, Pop2 = Pop2_new, Fst) %>%
  distinct(Pop1, Pop2, .keep_all = TRUE) %>%  # Remove duplicates
  mutate(
    Pop1 = factor(Pop1, levels = altitude_order),
    Pop2 = factor(Pop2, levels = altitude_order)
  ) %>%
  filter(as.numeric(Pop1) < as.numeric(Pop2))  # Keep only true lower triangle

# Plot with axes ordered low→high and dropped empty levels
p = ggplot(Fst_lower, aes(x = Pop2, y = Pop1)) +
  	geom_tile(aes(fill = Fst), color = "white") +
  	geom_text(aes(label = round(Fst, 3)), color = "white", size = 4) +
  	scale_fill_gradient(low = "#99D8C9", high = "#00441B", "Fst") +
  	scale_x_discrete(drop = TRUE, limits = altitude_order) +  # Hide empty columns
  	scale_y_discrete(drop = TRUE, limits = altitude_order) +  # Hide empty rows
  	labs(title = "Fst between Isophya populations") +
  	theme_minimal() +
  	theme(
    	axis.text.x = element_text(size = 12, colour = "black", face = "bold"),
    	axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
    	legend.title = element_text(size = 11, face = "bold"),  # Larger legend title
    	legend.text = element_text(size = 12),  # Larger legend labels
    	panel.grid = element_blank(),
    	axis.title = element_blank(),
    	plot.title = element_text(hjust = 0.5)
  	)

ggsave(
  "isophya71.Fst.pdf",
  p,
  width  = 250,
  height = 210,
  units  = "mm",
  dpi    = 300
)
