import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('markers') 
args = parser.parse_args()

myObject =  args.myObject
markers = args.markers
print(markers)

parts = myObject.split("_")
plot_name = parts[0] + ".png"


with open(markers, 'r') as file:
    marker_genes = [line.strip() for line in file.readlines()]

print(marker_genes)

combined_adata = sc.read_h5ad(myObject, backed="r")


sc.pl.dotplot(combined_adata, 
              var_names=marker_genes, 
              groupby='celltype',  
              save=plot_name) 

sc.pl.violin(
    combined_adata, 
    keys=marker_genes,
    groupby='celltype', rotation =90, 
    save=plot_name, stripplot=True,  
    show=True 
)



for gene in marker_genes:
    sc.pl.scatter(
        combined_adata,
        color=gene,  basis='umap', 
        save=plot_name + f"_{gene}_featureplot.png"
    )


