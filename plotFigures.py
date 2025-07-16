import scanpy as sc
import sys
import importlib_metadata
import anndata
import argparse
import scanpy.external as sce

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser(description="Plot expression for microglia genes.")
parser.add_argument("h5ad_file", help="Path to the .h5ad file")
args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad_file)

# Plot expression for il10, il10ra, il10rb
sc.pl.umap(adata, color=['il10', 'il10ra', 'il10rb'], save='_IL10_genes.png',use_raw=True) 

genes_to_plot = ['tnfa', 'tnfb', 'il1b']

# UMAP plots showing expression of each gene
sc.pl.umap(adata, color=genes_to_plot, 
           title=[f'Expression of {g}' for g in genes_to_plot], 
           save='_TNF_IL1b_expression.png',use_raw=True)


# Plot violin grouped by cell type
sc.pl.violin(
    adata,
    keys=genes_to_plot,
    groupby="celltype",  
    rotation = 90, 
    save="_TNF_IL1B_violin.png",use_raw=True)

adata_filtered = adata[adata.obs["renamed_samples"].isin(["Control", "LightDamaged"])].copy()

sc.pp.neighbors(adata_filtered)
sc.tl.umap(adata_filtered)


sc.pl.umap(
    adata_filtered,
    color="celltype",
    title="UMAP - Retinal Cell Clusters (Control & Light-Damaged)",
    frameon=False,
    legend_loc='on data',
    save="_controlandlightdamagedON.png"
)

sc.pl.umap(
    adata_filtered,
    color="celltype",
    title="UMAP - Retinal Cell Clusters (Control & Light-Damaged)",
    frameon=False,
    save="_controlandlightdamaged.png"
)
