import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np



sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('annotations')
args = parser.parse_args()

myObject =  args.myObject
annot_file = args.annotations

parts = myObject.split("_")
newObject = "annotated_" + parts[1]

fname = parts[1].split(".")
sample = fname[0]
combined_adata = sc.read_h5ad(myObject, backed="r")


cluster_to_celltype_dict = {}
with open(annot_file, "r") as f:
    for line in f:
        cluster, celltype = line.strip().split()
        cluster_to_celltype_dict[cluster] = celltype

cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}

combined_adata.obs["celltype"] = combined_adata.obs["leiden"].map(cluster_to_celltype_dict)


figure_name = sample +"_annotationsON.png"
combined_adata.obs_names_make_unique()
sc.pl.umap(combined_adata, color='celltype',legend_loc="on data", save=figure_name)
figure_name = sample +"_annotations.png"
sc.pl.umap(combined_adata, color='celltype',save=figure_name)

combined_adata.write(newObject, compression="gzip")
