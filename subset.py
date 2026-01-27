import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np
import os

# Fix for importlib.metadata
sys.modules['importlib.metadata'] = importlib_metadata

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('celltype')
args = parser.parse_args()
myObject = args.myObject
celltype = args.celltype

# Construct new object name
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = f"{celltype}_{base_name}.h5ad"

# Read the original object
adata = sc.read_h5ad(myObject)

# Subset the adata for the celltype of interest
subset = adata[adata.obs['celltype'].isin([celltype]), :].copy()

# Compute PCA, neighbors, UMAP, and Leiden clustering
sc.pp.pca(subset)
sc.pp.neighbors(subset)
sc.tl.umap(subset)
sc.tl.leiden(subset, resolution=1.0)

# Ensure output directory exists
os.makedirs("figures", exist_ok=True)

# UMAP plot for subset colored by Leiden
plot_name = f"umap{celltype}.png"
fig = sc.pl.umap(subset, color='leiden', legend_loc='on data', show=False, return_fig=True)
fig.set_size_inches(12, 12)
plt.savefig(f"figures/{plot_name}", dpi=600, bbox_inches="tight")
plt.close(fig)

# Ensure 'renamed_samples' is categorical in original adata
adata.obs['renamed_samples'] = adata.obs['renamed_samples'].astype('category')
renamed_samples = adata.obs['renamed_samples'].cat.categories

# Generate UMAPs for each sample individually
for sample in renamed_samples:
    sample_adata = adata[adata.obs['renamed_samples'] == sample].copy()
    sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].astype('category')
    sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].cat.remove_unused_categories()
    fig = sc.pl.umap(
        sample_adata,
        color='renamed_samples',
        title=f"Sample: {sample}",
        size=20,
        show=False,
        return_fig=True
    )
    fig.set_size_inches(12, 12)
    fig.savefig(f"figures/umap_{celltype}_{sample}.png", dpi=600, bbox_inches='tight')
    plt.close(fig)

# Plot all samples together
fig = sc.pl.umap(
    adata,
    color='renamed_samples',
    size=2,
    show=False,
    return_fig=True
)
fig.set_size_inches(12, 12)
fig.savefig(f"figures/umap_merged_{celltype}.png", dpi=600, bbox_inches='tight')
plt.close(fig)

# Finalize subset for saving: remove unused categories and raw reference
for col in subset.obs.select_dtypes(include='category').columns:
    subset.obs[col] = subset.obs[col].cat.remove_unused_categories()
if subset.raw is not None:
    subset.raw = None

# If file exists, remove it (fix for older Scanpy/AnnData versions)
if os.path.exists(newObject):
    os.remove(newObject)

# Save the subset
subset.write(newObject, compression="gzip")
print(f"Subset saved as {newObject}")

