#!/bin/bash

# Activate the r-heatmap environment
echo "Activating r-heatmap environment..."
micromamba activate r-heatmap

# Run the R script
echo "Running heatmap generation script..."
cd "$(dirname "$0")/scripts"
Rscript heatmap_final.R

echo ""
echo "✓ Complete!"
echo ""
echo "Output files:"
echo "  - PNG: results/heatmap_top10_genes.png"
echo "  - PDF: results/heatmap_top10_genes.pdf"
