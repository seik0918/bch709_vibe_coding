#!/usr/bin/env Rscript

# Install packages if needed
packages <- c("ggplot2", "dplyr", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# 1. Read the data
data <- read.delim("../data/gasch2000.txt", row.names = 1)

# 2. Check the structure
cat("Data dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
cat("First few rows and columns:\n")
print(head(data[,1:5]))

# 3. Check data range to determine if already log-scale
cat("\nData range (first 100 genes):\n")
cat("Min:", min(data[1:100,], na.rm=TRUE), "\n")
cat("Max:", max(data[1:100,], na.rm=TRUE), "\n")
cat("Values contain negatives (indicator of log-scale): ", any(data < 0, na.rm=TRUE), "\n")

# Since data contains negative values, it's already in log-scale (log fold-change)
# Do NOT apply additional log transformation

# 4. Extract gene expression values (exclude UID, NAME, GWEIGHT columns)
# The file has 3 annotation columns: UID (row names), NAME, GWEIGHT
# We keep all numeric expression columns starting from column 3 (GWEIGHT removed in practice)

# Check if column names include "NAME" and "GWEIGHT"
cat("\nFirst column names:", head(colnames(data)), "\n")

# Extract numeric expression data (skip NAME and GWEIGHT which are columns 1-2)
expr_data <- data[, 3:ncol(data)]
cat("Expression data dimensions:", nrow(expr_data), "x", ncol(expr_data), "\n")

# 5. Compute CV (coefficient of variation) per gene
# CV = sd / |mean| (using absolute mean to avoid division issues)
compute_cv <- function(x) {
  x_clean <- x[!is.na(x)]
  if(length(x_clean) == 0) return(NA)
  
  mean_val <- mean(x_clean)
  sd_val <- sd(x_clean)
  
  # Use absolute mean to handle cases where mean is near zero
  if(abs(mean_val) < 1e-10) return(NA)
  
  return(sd_val / abs(mean_val))
}

cv_values <- apply(expr_data, 1, compute_cv)

# Remove genes with NA CV
cv_values_clean <- cv_values[!is.na(cv_values)]
cat("Genes with valid CV:", length(cv_values_clean), "\n")

# 6. Get top 10 genes by CV
top10_idx <- names(sort(cv_values_clean, decreasing = TRUE)[1:10])
cat("\nTop 10 genes by coefficient of variation:\n")
print(sort(cv_values_clean, decreasing = TRUE)[1:10])

# Extract top 10 gene expression data
top10_data <- expr_data[top10_idx, ]

# 7. Reshape to long format for ggplot
library(dplyr)
library(tidyr)
top10_long <- top10_data %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "condition", values_to = "expression")

# 8. Create the heatmap
heatmap_plot <- ggplot(data = top10_long, aes(x = gene, y = condition, fill = expression)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_distiller(palette = "PuGn") +
  labs(
    title = "gene top 10",
    x = "gene",
    y = "conditions",
    fill = "Log Expression"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = "black", size = 1),
    plot.title = element_text(family = "Times New Roman", size = 11, hjust = 0.5),
    axis.title.x = element_text(family = "Times New Roman", size = 11),
    axis.title.y = element_text(family = "Times New Roman", size = 11),
    axis.text.x = element_text(family = "Times New Roman", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 9),
    legend.position = "top",
    legend.title = element_text(family = "Times New Roman", size = 11),
    legend.text = element_text(family = "Times New Roman", size = 10),
    panel.grid = element_blank()
  )

# 9. Export as high-resolution image (PNG format)
output_file <- "../results/heatmap_top10_genes.png"
ggsave(
  filename = output_file,
  plot = heatmap_plot,
  width = 5,
  height = 5,
  dpi = 300,
  units = "in",
  bg = "white"
)

cat("\nHeatmap saved to:", output_file, "\n")

# Also save as PDF for vector format
pdf_file <- "../results/heatmap_top10_genes.pdf"
ggsave(
  filename = pdf_file,
  plot = heatmap_plot,
  width = 5,
  height = 5,
  units = "in",
  bg = "white"
)

cat("Heatmap also saved to:", pdf_file, "\n")
