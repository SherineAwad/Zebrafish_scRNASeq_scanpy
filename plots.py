import scanpy as sc
import sys
import importlib_metadata
import anndata 
import argparse
import scanpy.external as sce
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os 

sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

adata = sc.read(myObject)
adata.var.head()
adata.obs.head()
adata.obs['leiden'].head()


figurename1 = "figures/sample_"+base_name+"_qc_violin.png"
sc.pl.violin(
    adata,
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='renamed_samples',     # change if you use another cluster label
    jitter=0.4,
    rotation=45,
    multi_panel=True,
    show=False
)
plt.savefig(figurename1, dpi=300, bbox_inches="tight")
plt.close()


figurename2 = "figures/celltype_"+base_name+"_qc_violin.png"
sc.pl.violin(
    adata,
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='celltype',     # change if you use another cluster label
    jitter=0.4,
    rotation=45,
    multi_panel=True,
    show=False
)
plt.savefig(figurename2, dpi=300, bbox_inches="tight")
plt.close()


# Validate required columns
if 'renamed_samples' not in adata.obs or 'celltype' not in adata.obs:
    raise ValueError("Your .h5ad file must contain 'renamed_samples' and 'celltype' in .obs")

# Sanitize values just in case
adata.obs['renamed_samples'] = adata.obs['renamed_samples'].astype(str).str.strip()
adata.obs['celltype'] = adata.obs['celltype'].astype(str).str.strip()


celltype_order = ['MG', 'MGPC', 'PR precursors', 'Rod', 'Cones', 'BC', 'AC', 'HC', 'RGC', 'Microlia','Microglia_ImmuneCells','RPE', 'Melanocyte','Endothelial','Perycites','Oligodenrocyte']


# Create dataframe of counts
df = (
    adata.obs[['renamed_samples', 'celltype']]
    .value_counts()
    .reset_index(name='count')
)

# Normalize to get proportions *per sample*
df['fraction'] = df['count'] / df.groupby('renamed_samples')['count'].transform('sum')

# Pivot for plotting: rows = samples, columns = celltypes
pivot_df = df.pivot(index='renamed_samples', columns='celltype', values='fraction').fillna(0)

# Ensure columns are in the desired order (filter out missing ones just in case)
ordered_columns = [ct for ct in celltype_order if ct in pivot_df.columns]
pivot_df = pivot_df[ordered_columns]

# Extract color mapping from adata if available
if 'celltype_colors' in adata.uns:
    # Ensure celltype is categorical with matching order
    if pd.api.types.is_categorical_dtype(adata.obs['celltype']):
        categories = adata.obs['celltype'].cat.categories
    else:
        categories = sorted(pivot_df.columns.tolist())  # fallback

    # Build color mapping dictionary
    celltype_colors = dict(zip(categories, adata.uns['celltype_colors']))

    # Reorder the color list to match new ordered_columns
    colors = [celltype_colors[ct] for ct in ordered_columns]
else:
    raise ValueError("'celltype_colors' not found in adata.uns. Please ensure UMAP plot colors are set.")

# Plot
ax = pivot_df.plot(kind='bar', stacked=True, figsize=(12, 6), color=colors)
plt.ylabel("Fraction of cells")
plt.xlabel("Sample")
plt.title("Cell Type Contribution per Sample")
plt.xticks(rotation=45, ha='right')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save figure
plt.savefig("figures/Restacked_bar_sample_by_celltype.png", dpi=300)
plt.close()


