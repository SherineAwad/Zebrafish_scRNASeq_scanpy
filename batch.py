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
newObject = "corrected_" + myObject

combined_adata = sc.read(myObject)



sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)

# HVG selection first

sc.pp.highly_variable_genes(combined_adata)
# Now save .raw — it includes normalized + log1p + HVG flags
combined_adata.raw = combined_adata.copy()

# Continue as usual
sc.pp.scale(combined_adata)
sc.tl.pca(combined_adata, svd_solver='arpack')
sce.pp.harmony_integrate(combined_adata, key='renamed_samples')
combined_adata.obsm['X_pca'] = combined_adata.obsm['X_pca_harmony']
sc.pp.neighbors(combined_adata, random_state=0)
sc.tl.umap(combined_adata, random_state=0)
sc.tl.leiden(combined_adata, resolution=1.0, random_state=0)



sc.pl.umap(combined_adata, color='renamed_samples', size=2, save='_Harmonyzebrafishes.png')

samples = combined_adata.obs['renamed_samples'].unique()

for sample in samples:
    sc.pl.umap(
        combined_adata[combined_adata.obs['renamed_samples'] == sample],
        color='renamed_samples',
        title=f"Sample: {sample}",
        size=20,
        save=f"_reclustered_{sample}.png",
        show=False
    )


combined_adata.write(newObject)



