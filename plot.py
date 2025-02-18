import scanpy as sc
import sys
import importlib_metadata
import anndata 
import argparse
import scanpy.external as sce

sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

combined_adata = sc.read(myObject)

combined_adata.var.head()
combined_adata.obs.head()
combined_adata.obs['leiden'].head()

sc.pl.violin(
    combined_adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="QCPerCluster.png", rotation=90)

sc.pl.violin(
    combined_adata,
    ["pct_counts_mt"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="mtPerCluster.png",rotation=90)

sc.pl.violin(
    combined_adata,
    ["total_counts"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="totalCountsPerCluster.png",rotation=90)


sc.pl.violin(
    combined_adata,
    ["n_genes_by_counts"],  # Cell-level attributes to plot
    jitter=0.4,
    multi_panel=True,
    groupby="leiden",save="geneCountsPerCluster.png",rotation=90)








