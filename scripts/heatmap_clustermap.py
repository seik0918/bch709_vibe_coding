#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

print("Loading libraries...")

# 1. Read the data
print("Reading data from ./data/gasch2000.txt")
df = pd.read_csv('./data/gasch2000.txt', sep='\t', index_col=0)

print(f"Data shape: {df.shape}")
print(f"First few columns: {df.columns[:5].tolist()}")
print(f"First few rows:")
print(df.iloc[:3, :3])

# 2. Extract numeric expression columns (skip NAME, GWEIGHT columns)
expr_data = df.iloc[:, 2:]
print(f"\nExpression data shape: {expr_data.shape}")

# 3. Compute CV (coefficient of variation) per gene
print("\nComputing CV for all genes...")
cv_values = expr_data.std(axis=1) / np.abs(expr_data.mean(axis=1))
cv_values = cv_values.replace([np.inf, -np.inf], np.nan)
cv_values = cv_values.dropna()

print(f"Genes with valid CV: {len(cv_values)}")

# Get top 10
top10_genes = cv_values.nlargest(10).index.tolist()

print("\nTop 10 genes by CV:")
for i, (gene, cv) in enumerate(cv_values.nlargest(10).items(), 1):
    print(f"{i}. {gene}: {cv:.4f}")

# 4. Extract top 10 genes, first 30 conditions
print("\nExtracting top 10 genes × first 30 conditions...")
top10_data = expr_data.loc[top10_genes, :].iloc[:, :30]

print(f"Final data shape: {top10_data.shape}")

# 5. Clean up index names (remove extra info)
top10_data.index = [name.split('..')[0] if '..' in name else name for name in top10_data.index]

# 6. Create the clustermap
print("\nCreating clustermap...")
g = sns.clustermap(
    top10_data,
    cmap='RdBu_r',              # Red-Blue diverging palette
    center=0,                   # Force white to be at 0
    annot=False,                # Keep clean
    figsize=(14, 8),            # Size for visibility
    dendrogram_ratio=0.15,      # Small dendrograms
    cbar_pos=(0.02, 0.8, 0.03, 0.15),  # Legend position
    linewidths=0.5,             # Add gridlines
    linecolor='gray',
    method='average',           # Linkage method for clustering
    metric='euclidean'          # Distance metric
)

# 7. Style improvements
g.ax_heatmap.set_xlabel('Conditions', fontsize=12, fontweight='bold')
g.ax_heatmap.set_ylabel('Genes', fontsize=12, fontweight='bold')
g.fig.suptitle('Top 10 Genes Heatmap (First 30 Conditions)', 
               fontsize=14, fontweight='bold', y=0.98)

# Rotate labels
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=11)

# Adjust layout
g.fig.tight_layout()

# 8. Create results directory if needed
import os
os.makedirs('./results', exist_ok=True)

# 9. Save outputs
print("\nSaving output files...")

png_file = './results/heatmap_clustermap_top10.png'
g.fig.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved PNG: {png_file}")

pdf_file = './results/heatmap_clustermap_top10.pdf'
g.fig.savefig(pdf_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved PDF: {pdf_file}")

plt.close()

print("\n✓ Done!")
