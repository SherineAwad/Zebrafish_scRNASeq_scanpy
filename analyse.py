import scanpy as sc
import sys
import importlib_metadata

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")  
newObject = "analysed_" + parts[1]

combined_adata = sc.read(myObject)

sc.pp.highly_variable_genes(combined_adata, flavor='seurat', n_top_genes=2000)

sc.pp.scale(combined_adata, max_value=10)

sc.tl.pca(combined_adata, svd_solver='arpack')

sc.pp.neighbors(combined_adata)

sc.tl.umap(combined_adata)

sc.pl.umap(combined_adata, color='sample', size=2, save='_zebrafishes.png')

combined_adata.write(newObject) 
