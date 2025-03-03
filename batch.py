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
parts = myObject.split("_")
newObject = "corrected_" + parts[1]

combined_adata = sc.read(myObject)


sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)
sc.pp.highly_variable_genes(combined_adata)
sc.pp.scale(combined_adata)

sc.tl.pca(combined_adata, svd_solver='arpack')
sce.pp.harmony_integrate(combined_adata, key='renamed_samples')

combined_adata.obsm['X_pca'] = combined_adata.obsm['X_pca_harmony']
sc.pp.neighbors(combined_adata, random_state=0)
sc.tl.umap(combined_adata)

sc.pl.umap(combined_adata, color='renamed_samples', size=2, save='_Harmonyzebrafishes.png')

combined_adata.write(newObject)
