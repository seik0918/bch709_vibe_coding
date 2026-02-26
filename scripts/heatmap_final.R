#!/usr/bin/env Rscript

# Force binary packages only, no source compilation
options(pkgType = "binary")

# Try to install ggplot2 with binary only
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "http://cran.r-project.org", type = "binary")
  library("ggplot2")
}

cat("ggplot2 loaded successfully\n")

# 1. Read the data
cat("Reading data from ./data/gasch2000.txt\n")
data <- read.delim("./data/gasch2000.txt", row.names = 1)

cat("Data dimensions:", nrow(data), "rows x", ncol(data), "columns\n")

# 2. Extract expression values (skip annotation columns 1-2)
expr_data <- data[, 3:ncol(data)]
cat("Expression data dimensions:", nrow(expr_data), "x", ncol(expr_data), "\n")

# 3. Compute CV (coefficient of variation)
compute_cv <- function(x) {
  x_clean <- x[!is.na(x)]
  if(length(x_clean) < 2) return(NA)
  
  m <- mean(x_clean)
  s <- sd(x_clean)
  if(abs(m) < 1e-10) return(NA)
  
  return(s / abs(m))
}

cat("Computing CV for all genes...\n")
cv_values <- apply(expr_data, 1, compute_cv)
cv_valid <- cv_values[!is.na(cv_values)]

cat(length(cv_valid), "genes with valid CV\n")

# Get top 10
top10_idx <- names(sort(cv_valid, decreasing = TRUE)[1:10])

cat("\nTop 10 genes by CV:\n")
tmp <- sort(cv_valid, decreasing = TRUE)[1:10]
for (i in seq_along(tmp)) {
  cat(sprintf("%d. %s: %.4f\n", i, names(tmp)[i], tmp[i]))
}

# Extract top 10 data
top10_data <- expr_data[top10_idx, ]

# Keep only first 30 columns (conditions)
cat("Using first 30 columns/conditions...\n")
top10_data <- top10_data[, 1:min(30, ncol(top10_data))]

# Load pheatmap for better visualization (should be available from conda)
if (!require("pheatmap", quietly = TRUE)) {
  cat("pheatmap package not found; attempting install from CRAN\n")
  install.packages("pheatmap", repos = "http://cran.r-project.org")
  library(pheatmap)
}

# Scale rows (z-score) to improve contrast across genes
scaled_data <- t(scale(t(top10_data)))
scaled_data[is.na(scaled_data)] <- 0

# Determine clustering order if needed
row_clust <- hclust(dist(scaled_data))
col_clust <- hclust(dist(t(scaled_data)))

# Reorder matrix for ggplot version
ordered_genes <- row_clust$labels[row_clust$order]
ordered_cols <- col_clust$labels[col_clust$order]

ordered_matrix <- scaled_data[ordered_genes, ordered_cols]

# Reshape to long format for ggplot
cat("\nReshaping data to long format...\n")
long_data <- data.frame()
for (i in 1:nrow(ordered_matrix)) {
  for (j in 1:ncol(ordered_matrix)) {
    val <- ordered_matrix[i, j]
    if (!is.na(val)) {
      long_data <- rbind(long_data, data.frame(
        gene = rownames(ordered_matrix)[i],
        condition = colnames(ordered_matrix)[j],
        expression = val,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("Long format data:", nrow(long_data), "rows\n")

# Create heatmap
cat("\nCreating heatmap...\n")

heatmap_plot <- ggplot(long_data, aes(x = gene, y = condition, fill = expression)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  ggtitle("Top 10 genes (scaled & clustered)") +
  xlab("gene") +
  ylab("conditions") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

# Save the ggplot version first
gg_png <- "./results/heatmap_top10_genes_first30cols.png"
cat("\nSaving ggplot PNG...\n")
ggsave(gg_png, heatmap_plot, width = 5, height = 5, dpi = 300, units = "in", bg = "white")
cat("Saved to:", gg_png, "\n")

# also save ggplot PDF
gg_pdf <- "./results/heatmap_top10_genes_first30cols.pdf"
cat("Saving ggplot PDF...\n")
ggsave(gg_pdf, heatmap_plot, width = 5, height = 5, units = "in", bg = "white")
cat("Saved to:", gg_pdf, "\n")

# Now create a nicer clustered heatmap with pheatmap
# transpose matrix so genes appear on x-axis and conditions on y-axis
cat("\nGenerating clustered heatmap with pheatmap (genes on x-axis)...\n")
png_clust <- "./results/heatmap_clustered_top10.png"
pdf_clust <- "./results/heatmap_clustered_top10.pdf"

# we'll use t(scaled_data) for orientation flip
flipped <- t(scaled_data)

pheatmap(flipped,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         border_color = "black",
         fontsize_row = 10,
         fontsize_col = 8,
         main = "Top 10 genes (scaled & clustered)",
         filename = png_clust,
         width = 5,
         height = 5)

cat("Saved clustered PNG to:", png_clust, "\n")

# PDF version
pheatmap(flipped,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         border_color = "black",
         fontsize_row = 10,
         fontsize_col = 8,
         main = "Top 10 genes (scaled & clustered)",
         filename = pdf_clust,
         width = 5,
         height = 5)
cat("Saved clustered PDF to:", pdf_clust, "\n")

# Create results directory if it doesn't exist
if (!dir.exists("./results")) {
  dir.create("./results", recursive = TRUE)
  cat("Created ./results directory\n")
}

# Save PNG
cat("\nSaving PNG...\n")
png_file <- "./results/heatmap_top10_genes_first30cols.png"
ggsave(png_file, heatmap_plot, width = 5, height = 5, dpi = 300, units = "in", bg = "white")
cat("Saved to:", png_file, "\n")

# Save PDF
cat("Saving PDF...\n")
pdf_file <- "./results/heatmap_top10_genes_first30cols.pdf"
ggsave(pdf_file, heatmap_plot, width = 5, height = 5, units = "in", bg = "white")
cat("Saved to:", pdf_file, "\n")

cat("\n✓ Done!\n")
