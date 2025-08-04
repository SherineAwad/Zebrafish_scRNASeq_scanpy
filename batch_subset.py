import scanpy as sc
import sys
import importlib_metadata
import anndata
import argparse
import scanpy.external as sce
import os 
import numpy as np 

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
sce.pp.harmony_integrate(combined_adata, key='renamed_samples')
combined_adata.obsm['X_pca'] = combined_adata.obsm['X_pca_harmony']
sc.pp.neighbors(combined_adata, random_state=0)
sc.tl.umap(combined_adata, random_state=0)
sc.tl.leiden(combined_adata, resolution=1.0, random_state=0)



sc.pl.umap(combined_adata, color='renamed_samples', size=2, save=f'_{base_name}_Harmony.png')

samples = combined_adata.obs['renamed_samples'].unique()

for sample in samples:
    sc.pl.umap(
        combined_adata[combined_adata.obs['renamed_samples'] == sample],
        color='renamed_samples',
        title=f"Sample: {sample}",
        size=20,
        save=f"_{base_name}_Harmony_{sample}.png",
        show=False
    )


combined_adata.write(newObject)



