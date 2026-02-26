#!/usr/bin/env Rscript

# Install packages with binary versions
packages <- c("ggplot2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org", type = "binary")
    library(pkg, character.only = TRUE)
  }
}

# 1. Read the data
data <- read.delim("../data/gasch2000.txt", row.names = 1)

cat("Data dimensions:", nrow(data), "rows x", ncol(data), "columns\n")
cat("First few column names:", head(colnames(data)), "\n")

# 2. Extract expression values (skip NAME and GWEIGHT columns - columns 1 and 2)
expr_data <- data[, 3:ncol(data)]
cat("Expression data dimensions after removing annotations:", nrow(expr_data), "x", ncol(expr_data), "\n")

# 3. Compute CV (coefficient of variation) per gene
compute_cv <- function(x) {
  x_clean <- x[!is.na(x)]
  if(length(x_clean) == 0) return(NA)
  
  mean_val <- mean(x_clean)
  sd_val <- sd(x_clean)
  
  if(abs(mean_val) < 1e-10) return(NA)
  
  return(sd_val / abs(mean_val))
}

cv_values <- apply(expr_data, 1, compute_cv)
cv_values_clean <- cv_values[!is.na(cv_values)]

cat("Genes with valid CV:", length(cv_values_clean), "\n")

# 4. Get top 10 genes by CV
top10_names <- names(sort(cv_values_clean, decreasing = TRUE)[1:10])

cat("\nTop 10 genes by coefficient of variation:\n")
print(data.frame(
  gene = top10_names,
  cv = sort(cv_values_clean, decreasing = TRUE)[1:10]
))

# Extract top 10 data
top10_data <- expr_data[top10_names, ]

# 5. Reshape to long format using base R (no dplyr needed)
long_data <- data.frame()
for (i in 1:nrow(top10_data)) {
  gene_name <- rownames(top10_data)[i]
  for (j in 1:ncol(top10_data)) {
    cond_name <- colnames(top10_data)[j]
    expr_value <- top10_data[i, j]
    long_data <- rbind(long_data, data.frame(
      gene = gene_name,
      condition = cond_name,
      expression = expr_value
    ))
  }
}

# Remove rows with NA expression values
long_data <- long_data[!is.na(long_data$expression), ]

cat("\nLong format data dimensions:", nrow(long_data), "rows\n")

# 6. Create the heatmap
heatmap_plot <- ggplot(data = long_data, aes(x = gene, y = condition, fill = expression)) +
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
    plot.title = element_text(family = "Times New Roman", size = 11, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(family = "Times New Roman", size = 11),
    axis.title.y = element_text(family = "Times New Roman", size = 11),
    axis.text.x = element_text(family = "Times New Roman", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(family = "Times New Roman", size = 9),
    legend.position = "top",
    legend.title = element_text(family = "Times New Roman", size = 11),
    legend.text = element_text(family = "Times New Roman", size = 10),
    panel.grid = element_blank()
  )

# 7. Export as high-resolution PNG
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

# Also save as PDF
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
cat("\nDone!\n")
