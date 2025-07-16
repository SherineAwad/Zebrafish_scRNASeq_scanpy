import scanpy as sc
import sys
import importlib_metadata
import anndata 
import argparse
import scanpy.external as sce
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

adata = sc.read(myObject)
'''
adata.var.head()
adata.obs.head()
adata.obs['leiden'].head()

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="QCPerCluster.png", rotation=90)

sc.pl.violin(
    adata,
    ["pct_counts_mt"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="mtPerCluster.png",rotation=90)

sc.pl.violin(
    adata,
    ["total_counts"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="totalCountsPerCluster.png",rotation=90)


sc.pl.violin(
    adata,
    ["n_genes_by_counts"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="geneCountsPerCluster.png",rotation=90)
'''


# Validate required columns
if 'renamed_samples' not in adata.obs or 'celltype' not in adata.obs:
    raise ValueError("Your .h5ad file must contain 'renamed_samples' and 'celltype' in .obs")

# Create dataframe of counts
df = (
    adata.obs[['celltype', 'renamed_samples']]
    .value_counts()
    .reset_index(name='count')
)

# Normalize to get proportions *per cell type*
df['fraction'] = df['count'] / df.groupby('celltype')['count'].transform('sum')

# Pivot for plotting
pivot_df = df.pivot(index='celltype', columns='renamed_samples', values='fraction').fillna(0)

# Plot
ax = pivot_df.plot(kind='bar', stacked=True, figsize=(12, 6), colormap='tab20')
plt.ylabel("Fraction of cells")
plt.xlabel("Cell Type")
plt.title("Sample contribution per Cell Type")
plt.xticks(rotation=45, ha='right')
plt.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save figure
plt.savefig("figures/Restacked_bar_celltype_by_sample.png", dpi=300)
plt.close()

print("✅ Saved: figures/stacked_bar_celltype_by_sample.png")

