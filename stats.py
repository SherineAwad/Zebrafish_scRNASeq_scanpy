import scanpy as sc
import sys
import importlib_metadata
import anndata
import argparse
import scanpy.external as sce
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
base_name = os.path.splitext(os.path.basename(myObject))[0]

adata = sc.read(myObject)

# Validate required columns
if 'renamed_samples' not in adata.obs or 'celltype' not in adata.obs:
    raise ValueError("Your .h5ad file must contain 'renamed_samples' and 'celltype' in .obs")

# Sanitize values
adata.obs['renamed_samples'] = adata.obs['renamed_samples'].astype(str).str.strip()
adata.obs['celltype'] = adata.obs['celltype'].astype(str).str.strip()

# Count cells per renamed_samples group
cell_counts = adata.obs['renamed_samples'].value_counts()
print("Number of cells per renamed_samples group:")
print(cell_counts)



marker_genes = {
"MG":	["aqp4"],
"MGPC":	["cdk1"],
"PR precursors": ["nr2e3"],
"Rod": ["guca1b"],
"Cones": ["gnat2"],
"BC":	["cabp5a"],
"AC": ["elavl3"],
"HC": 	["ompa"],
"RGC":	["isl2b"],
"Microglia_ImmuneCells": ["mpeg1.1"],
"RPE": 	["rpe65a"],
"Endothelial": ["tie1"],
"Pericytes": ["acta2"],
"Oligocytes": ["mbpa"],
"Melanocytes": ["mitfa"]
}
celltype_order = ['MG', 'MGPC', 'PR precursors', 'Rod', 'Cones', 'BC', 'AC', 'HC', 'RGC','Microglia_ImmuneCells','RPE', 'Melanocyte','Endothelial','Perycites','Oligodenrocyte']
adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'], categories=celltype_order, ordered=True)

figure_name = f"figures/dotplot_{base_name}_markerGenes.png"

'''
fig = sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="celltype",
    categories_order = celltype_order,
    standard_scale="var",
    figsize=(6,5),
    dot_max=1.0, dendrogram=False, show =False,return_fig=True)

fig.savefig(figure_name, dpi=600,bbox_inches="tight")

'''

figure_name = f"figures/dotplot_{base_name}_markerGenes.png"
fig = sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="celltype",
    categories_order=celltype_order,
    standard_scale="var",
    figsize=(6,5),
    dot_max=1.0,
    dendrogram=False,
    show=False,
    return_fig=True
)

fig.savefig(figure_name, dpi=600, bbox_inches="tight")


