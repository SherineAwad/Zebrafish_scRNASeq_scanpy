import decoupler as dc
import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import omnipath 


sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject

parts = myObject.split("_")
newObject = "enriched_" + parts[1]

fname = parts[1].split(".")
filename = fname[0] + "_pathways.png" 

combined_adata = sc.read_h5ad(myObject, backed="r")
combined_adata = combined_adata.to_memory()

dc.show_resources()

msigdb_zebrafish = dc.get_resource('MSigDB', organism='danio_rerio')

dc.run_ora(mat=combined_adata, net=msigdb_zebrafish, source='source', target='target')

combined_adata.obsm['pathway_activities'] = combined_adata.obsm['ora_estimate']

# Plot the UMAP with the pathway activities colored
sc.pl.umap(combined_adata, color=['pathway_activities'], cmap='coolwarm', save=filename)

# Save the modified AnnData object
combined_adata.write(newObject, compression="gzip")



