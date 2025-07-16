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

with open(markers, 'r') as f:
    marker_genes = [g.strip() for g in f if g.strip()]


present = [g for g in marker_genes if g in combined_adata.var_names]
missing = [g for g in marker_genes if g not in combined_adata.var_names]  
if missing:
    print(f"Warning: these genes were not found in the dataset and will be skipped:\n  {missing}")

if present:
    sc.pl.dotplot(combined_adata, var_names=present, groupby='celltype', save=plot_name)
    sc.pl.violin(combined_adata, keys=present, groupby='celltype',
                 rotation=90, stripplot=True, save=plot_name, show=True)
else:
    print("No marker genes found—skipping dotplot/violin.")

for gene in present:
    sc.pl.scatter(combined_adata,
                  color=gene,
                  basis='umap',
                  save=f"{plot_name}_{gene}_featureplot.png")





