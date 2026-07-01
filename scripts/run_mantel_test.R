# Load necessary library
library(tidyr)
library(dplyr)
library(vegan)

# Read the data
data <- read.table("/home/ismail/Research/isophya/isophya71.fst.dist.tsv", header = TRUE, sep = " ", stringsAsFactors = FALSE)
data$Linearized_Fst <- data$Fst / (1 - data$Fst)


# Create a function to generate a full matrix
generate_matrix <- function(data, value_col) {
  # Pivot data to wide format
  matrix <- data %>%
    select(Pop1, Pop2, all_of(value_col)) %>%
    pivot_wider(names_from = Pop2, values_from = all_of(value_col))
  
  # Convert the "Pop1" column to row names
  row.names(matrix) <- matrix$Pop1
  matrix <- matrix[, -1]  # Remove the "Pop1" column from the data frame
  
  # Fill missing values by mirroring the upper triangle to the lower triangle
  matrix <- as.matrix(matrix)
  matrix[lower.tri(matrix)] <- t(matrix)[lower.tri(matrix)]
  
  return(matrix)
}

# Generate the distance matrix
distance_matrix <- generate_matrix(data, "Distance")

# Generate the Fst matrix
linearized_fst_matrix <- generate_matrix(data, "Linearized_Fst")

# Save the matrices to files (optional)
write.table(distance_matrix, "/home/ismail/Research/isophya/distance_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(linearized_fst_matrix, "/home/ismail/Research/isophya/linearized_fst_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)


# Ensure matrices are symmetrical and of the same size
if (!all(dim(distance_matrix) == dim(linearized_fst_matrix))) {
  stop("The matrices must have the same dimensions.")
}

# Perform Mantel test using Pearson correlation
mantel_pearson <- mantel(distance_matrix, linearized_fst_matrix, method = "pearson", permutations = 999)

# Perform Mantel test using Spearman correlation
mantel_spearman <- mantel(distance_matrix, linearized_fst_matrix, method = "spearman", permutations = 999)

# Format the results into a data frame
mantel_results <- data.frame(
  Method = c("Pearson", "Spearman"),
  Statistic = c(mantel_pearson$statistic, mantel_spearman$statistic),
  P_value = c(mantel_pearson$signif, mantel_spearman$signif)
)

# Save the results to a TSV file
write.table(mantel_results, "/home/ismail/Research/isophya/mantel_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Print a message indicating success
print("Mantel test results saved to 'mantel_results.tsv'")
