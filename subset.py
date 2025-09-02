import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np
import os


sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('celltype') 
args = parser.parse_args()
myObject =  args.myObject

parts = myObject.split("_")  
celltype = args.celltype 
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject = f"{celltype}_{base_name}"
adata = sc.read_h5ad(myObject)
subset = adata[adata.obs['celltype'].isin([celltype]), :].copy()

sc.pp.pca(subset)
sc.pp.neighbors(subset)
sc.tl.umap(subset)
sc.tl.leiden(subset,resolution=1.0)

plot_name = "umap"+celltype+".png" 
fig = sc.pl.umap(subset, color='leiden', legend_loc='on data', show=False,return_fig=True) 
fig.set_size_inches(12, 12)
plt.savefig("figures/" + plot_name, dpi=600, bbox_inches="tight")
plt.close(fig)

# Ensure 'renamed_samples' is categorical
adata.obs['renamed_samples'] = adata.obs['renamed_samples'].astype('category')
renamed_samples = adata.obs['renamed_samples'].cat.categories
for sample in renamed_samples:
    sample_adata = adata[adata.obs['renamed_samples'] == sample].copy()
    sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].astype('category')
    sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].cat.remove_unused_categories()
    fig = sc.pl.umap(
        sample_adata,
        color='renamed_samples',
        title=f"Sample: {sample}",
        size=20,
        show=False, return_fig=True)
    fig.savefig(f"figures/umap_{celltype}_{sample}.png", dpi=600,bbox_inches='tight')
    fig.set_size_inches(12, 12)
plt.close(fig)

# Plot all samples together
fig = sc.pl.umap(
    adata,
    color='renamed_samples',
    size=2, show=False, return_fig=True) 
fig.set_size_inches(12, 12)
fig.savefig(f"figures/umap_merged_{celltype}.png", dpi=600,bbox_inches='tight')
plt.close(fig)

newObject = newObject +".h5ad"
subset.write(newObject, compression="gzip", overwrite=True)

