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

newObject = "annotated_" + myObject 

sample = "Zebrafishes" 
combined_adata = sc.read_h5ad(myObject)


cluster_to_celltype_dict = {}
with open(annot_file, "r") as f:
    for line in f:
        cluster, celltype = line.strip().split('\t')
        cluster_to_celltype_dict[cluster] = celltype

cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}

combined_adata.obs["celltype"] = combined_adata.obs["leiden"].map(cluster_to_celltype_dict)


figure_name = sample +"_annotationsON.png"
combined_adata.obs_names_make_unique()
sc.pl.umap(combined_adata, color='celltype',legend_loc="on data", save=figure_name)

figure_name = sample +"_annotations.png"
sc.pl.umap(combined_adata, color='celltype',save=figure_name)

unwanted_type = "Cones_MG_MGPC_PostMitotic"
combined_adata_filtered = combined_adata[combined_adata.obs["celltype"] != unwanted_type].copy()

# Plot UMAP with annotations
sc.pl.umap(combined_adata_filtered, color='celltype', legend_loc="on data", save=sample + "_NannotationsON.png")
sc.pl.umap(combined_adata_filtered, color='celltype', save=sample + "_Nannotations.png")


combined_adata_filtered.write(newObject, compression="gzip")
