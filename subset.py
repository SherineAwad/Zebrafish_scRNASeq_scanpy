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
parser.add_argument('celltype') 

args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")  

celltype = args.celltype 

combined_adata = sc.read_h5ad(myObject, backed="r")

subset = combined_adata[combined_adata.obs['celltype'].isin([celltype]), :]

subset_name = celltype +"_"+ parts[1] 

subset.write(subset_name, compression="gzip")

