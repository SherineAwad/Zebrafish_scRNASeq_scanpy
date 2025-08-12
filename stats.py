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


