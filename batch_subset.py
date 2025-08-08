import scanpy as sc
import sys
import importlib_metadata
import anndata
import argparse
import scanpy.external as sce
import os 
import numpy as np 
import matplotlib.pyplot as plt

sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
newObject = "corrected_" + myObject

base_name = os.path.splitext(os.path.basename(myObject))[0]

combined_adata = sc.read(myObject)

sc.pp.filter_cells(combined_adata, min_counts=1)

import scipy.sparse
if scipy.sparse.issparse(combined_adata.X):
    X_vals = combined_adata.X.data
else:
    X_vals = combined_adata.X

print("Before normalization:")
print("Min:", np.min(X_vals))
print("Max:", np.max(X_vals))
print("NaNs:", np.isnan(X_vals).any())


# Optional: Ensure dense matrix if needed
if scipy.sparse.issparse(combined_adata.X):
    combined_adata.X = combined_adata.X.toarray()

# Clean and normalize
combined_adata.X = np.clip(combined_adata.X, 0, None)
sc.pp.normalize_total(combined_adata)
sc.pp.log1p(combined_adata)
combined_adata.X = np.nan_to_num(combined_adata.X, nan=0.0, posinf=0.0, neginf=0.0)
# HVG detection and scaling
sc.pp.highly_variable_genes(combined_adata)
combined_adata.raw = combined_adata.copy()
sc.pp.scale(combined_adata)
sc.tl.pca(combined_adata, svd_solver='arpack')
# Perform batch correction using Harmony
sce.pp.harmony_integrate(
    combined_adata,
    key='renamed_samples',  
    # 🔸 key: str
    # This is the column in `combined_adata.obs` that identifies batches, 
    # such as sample names or experimental conditions.
    # Harmony will align the embeddings across the values in this column.

    theta=3,  
    # 🔸 theta: float (default = 2)
    # Controls the **strength of batch correction**.
    # Higher values (e.g., 5–10) make Harmony more aggressive in aligning batches,
    # possibly at the expense of biological variation.
    # Lower values (e.g., 0–1) make it more conservative — useful if you want to retain subtle differences.

    max_iter_harmony=20,  
    # 🔸 max_iter_harmony: int (default = 20)
    # Sets the **maximum number of optimization iterations** Harmony will run.
    # More iterations allow better convergence, especially on complex datasets.
    # Increase to 30–50 if Harmony stops before reaching convergence.

    epsilon_cluster=1e-6,  
    # 🔸 epsilon_cluster: float (default = 1e-5)
    # **Tolerance for convergence of the soft clustering step.**
    # Smaller values make convergence stricter (more accurate, slower).
    # Can usually be left at default unless tuning for performance.

    epsilon_harmony=1e-5,  
    # 🔸 epsilon_harmony: float (default = 1e-4)
    # **Tolerance for convergence of the batch correction step itself.**
    # Smaller values force Harmony to refine the correction more precisely.

    sigma=0.1,  
    # 🔸 sigma: float (default = 0.1)
    # Controls the **softness of cluster assignments** in Harmony's iterative updates.
    # Lower values (e.g., 0.05) make clustering harder (more distinct clusters).
    # Higher values (e.g., 0.2) soften boundaries between clusters.

    tau=0  
    # 🔸 tau: int (default = 0)
    # The number of **initial PCA components to exclude** from correction.
    # Use this if you know early PCs are dominated by technical noise or batch effect 
    # and want to avoid misleading Harmony in early updates.
    # Usually safe to keep as 0.
)

combined_adata.obsm['X_pca'] = combined_adata.obsm['X_pca_harmony']
sc.pp.neighbors(combined_adata, random_state=0)
sc.tl.umap(combined_adata, random_state=0)
sc.tl.leiden(combined_adata, resolution=1.0, random_state=0)



fig = sc.pl.umap(combined_adata, color='renamed_samples', size=2, show=False,return_fig=True)
fig.savefig(f'figures/{base_name}_Harmony.png', dpi=600,bbox_inches='tight')
plt.close(fig)

samples = combined_adata.obs['renamed_samples'].unique()

for sample in samples:
    fig = sc.pl.umap(
        combined_adata[combined_adata.obs['renamed_samples'] == sample],
        color='renamed_samples',
        title=f"Sample: {sample}",
        size=20,
        show=False, return_fig=True
    )
    fig.savefig(f"figures/{base_name}_Harmony_{sample}.png", dpi=600,bbox_inches='tight')
    plt.close(fig)
    

fig = sc.pl.umap(
    combined_adata,
    color='leiden',
    title="Leiden clusters",
    size=20,
    legend_loc="on data",
    show=False,
    return_fig=True
)
fig.savefig(f"figures/{base_name}_Harmony_leiden.png", dpi=600, bbox_inches='tight')
plt.close(fig) 

sc.pl.violin(
    combined_adata,
    keys=['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    groupby='leiden',
    jitter=0.4,
    rotation=45,
    multi_panel=False,
    show=False,
    save=f"_{base_name}_qc.png"
)


combined_adata.write(newObject)
