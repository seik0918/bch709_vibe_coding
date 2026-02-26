#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "gasch2000.txt"
output_path <- if (length(args) >= 2) args[[2]] else "results/gene_top10_heatmap.png"

download_url <- "https://www.shackett.org/files/gasch2000.txt"

resolve_input <- function(path) {
  if (file.exists(path)) return(path)
  alt <- file.path("data", basename(path))
  if (file.exists(alt)) return(alt)
  return(path)
}

input_path <- resolve_input(input_path)

if (!file.exists(input_path)) {
  message("Input not found. Downloading: ", download_url)
  status <- system2("curl", c("-L", "-o", input_path, download_url))
  if (!identical(status, 0L) || !file.exists(input_path)) {
    stop("Failed to download input file: ", input_path)
  }
}

message("Reading: ", input_path)
df <- read.delim(
  input_path,
  sep = "\t",
  header = TRUE,
  quote = "\"",
  comment.char = "#",
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

if (ncol(df) < 2) {
  stop("Input does not look like a matrix with gene IDs + expression columns.")
}

gene_ids <- trimws(as.character(df[[1]]))
gene_ids[gene_ids == "" | is.na(gene_ids)] <- paste0("gene_", which(gene_ids == "" | is.na(gene_ids)))
gene_ids <- make.unique(gene_ids)

# Keep only columns that are mostly numeric (expression columns), dropping annotation columns.
numeric_col <- vapply(df[-1], function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  mean(!is.na(x_num)) >= 0.8
}, logical(1))

expr <- as.data.frame(lapply(df[-1][, numeric_col, drop = FALSE], function(x) suppressWarnings(as.numeric(x))),
                      check.names = FALSE)
expr <- expr[, toupper(colnames(expr)) != "GWEIGHT", drop = FALSE]

if (ncol(expr) < 30) {
  stop("Fewer than 30 numeric condition columns found (", ncol(expr), ").")
}

rownames(expr) <- gene_ids
expr30 <- expr[, seq_len(30), drop = FALSE]
clean_conditions <- iconv(colnames(expr30), from = "", to = "ASCII//TRANSLIT", sub = "?")
clean_conditions[is.na(clean_conditions) | clean_conditions == ""] <- paste0(
  "condition_",
  which(is.na(clean_conditions) | clean_conditions == "")
)
clean_conditions <- make.unique(clean_conditions)
condition_labels <- sprintf("C%02d", seq_len(ncol(expr30)))
condition_map <- data.frame(label = condition_labels, condition = clean_conditions, stringsAsFactors = FALSE)
colnames(expr30) <- condition_labels

message("Structure check")
message("- Genes (rows): ", nrow(expr30))
message("- Numeric conditions used: ", ncol(expr30))
message("- Missing values: ", sum(is.na(expr30)))
message("- Contains negative values: ", any(expr30 < 0, na.rm = TRUE))
message("- Value range: [", paste(signif(range(expr30, na.rm = TRUE), 4), collapse = ", "), "]")
message("- CV filter: at least 24 non-NA values and abs(mean) >= 0.1")
message("- Condition key (C01..C30):")
print(condition_map)

# If negatives are present, data is likely already log-scale; otherwise apply log2(x + 1).
if (any(expr30 < 0, na.rm = TRUE)) {
  expr_log <- expr30
  message("Detected likely log-scale matrix. No additional log transform applied.")
} else {
  expr_log <- log2(expr30 + 1)
  message("No negatives detected. Applied log2(x + 1).")
}

compute_cv <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 24) return(NA_real_)
  m <- mean(x)
  if (abs(m) < 0.1) return(NA_real_)
  sd(x) / abs(m)
}

cv <- apply(expr_log, 1, compute_cv)
cv <- cv[!is.na(cv)]

if (length(cv) < 10) {
  stop("Fewer than 10 genes with valid CV after filtering.")
}

top10 <- names(sort(cv, decreasing = TRUE))[1:10]
message("Top 10 genes by CV:")
print(round(sort(cv, decreasing = TRUE)[1:10], 4))

top_mat <- expr_log[top10, , drop = FALSE]

# Base-R reshape to keep dependencies minimal.
plot_df <- data.frame(
  gene = rep(rownames(top_mat), times = ncol(top_mat)),
  condition = rep(colnames(top_mat), each = nrow(top_mat)),
  expression = as.vector(as.matrix(top_mat)),
  stringsAsFactors = FALSE
)

plot_df$gene <- factor(plot_df$gene, levels = top10)
plot_df$condition <- factor(plot_df$condition, levels = rev(colnames(top_mat)))

make_plot <- function(font_family) {
  ggplot(plot_df, aes(x = gene, y = condition, fill = expression)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = colorRampPalette(RColorBrewer::brewer.pal(9, "PuBuGn"))(256),
      na.value = "grey90"
    ) +
    labs(
      title = "gene top 10",
      x = "gene",
      y = "conditions",
      fill = "expression"
    ) +
    theme_classic(base_family = font_family, base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

out_dir <- dirname(output_path)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

device_fun <- NULL
if (requireNamespace("ragg", quietly = TRUE)) {
  device_fun <- ragg::agg_png
}

saved <- FALSE
for (font_try in c("Times New Roman", "Times", "sans")) {
  attempt <- try(
    ggsave(
      filename = output_path,
      plot = make_plot(font_try),
      width = 5,
      height = 5,
      dpi = 300,
      units = "in",
      bg = "white",
      device = device_fun
    ),
    silent = TRUE
  )
  if (!inherits(attempt, "try-error") && file.exists(output_path) && file.info(output_path)$size > 0) {
    saved <- TRUE
    message("Saved using font family: ", font_try)
    break
  }
}

if (!saved) {
  stop("Failed to save plot. In the conda environment, ensure r-ragg is installed.")
}
if (!file.exists(output_path)) {
  stop("Heatmap save failed unexpectedly.")
}

if (file.info(output_path)$size == 0) {
  stop("Heatmap file is empty.")
}

message("Saved heatmap: ", output_path)
