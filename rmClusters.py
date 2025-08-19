import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import numpy as np
import os

sys.modules['importlib.metadata'] = importlib_metadata

# ---------------------------
# Parse arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Remove specific clusters from an AnnData object and replot.")
parser.add_argument('myObject', help="Path to input h5ad file")
parser.add_argument('clusters_file', help="Path to text file with clusters to remove (one per line)")
args = parser.parse_args()

myObject = args.myObject
clusters_file = args.clusters_file

# ---------------------------
# Load AnnData
# ---------------------------
base_name = os.path.splitext(os.path.basename(myObject))[0]
newObject_name = "filtered_" + base_name
adata = sc.read_h5ad(myObject)

# ---------------------------
# Read clusters to remove from file
# ---------------------------
unwanted_clusters = []
with open(clusters_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line:  # ignore empty lines
            unwanted_clusters.append(int(line))

if len(unwanted_clusters) == 0:
    print("No clusters found in the file. Exiting.")
    sys.exit(1)

# ---------------------------
# Ensure 'leiden' column is int
# ---------------------------
adata.obs['leiden'] = adata.obs['leiden'].astype(int)

# ---------------------------
# Remove unwanted clusters using isin (same style as unwanted_types example)
# ---------------------------
adata_filtered = adata[~adata.obs["leiden"].isin(unwanted_clusters)].copy()
print(f"Removed clusters: {unwanted_clusters}")
print(f"Remaining clusters: {sorted(adata_filtered.obs['leiden'].unique())}")

# ---------------------------
# Plot
# ---------------------------
# Ensure leiden is categorical
adata_filtered.obs['leiden'] = adata_filtered.obs['leiden'].astype('category')
figure_name = "figures/"+ base_name +"_filteredannotationsON.png"
fig = sc.pl.umap(adata_filtered, color=["leiden"], save=figure_name, legend_loc="on data",show = False, return_fig=True)
fig.savefig(figure_name, dpi=600,bbox_inches='tight')
plt.close(fig)

# ---------------------------
# Save filtered object
# ---------------------------
adata_filtered.write_h5ad(newObject_name + ".h5ad")
print(f"Filtered AnnData saved as: {newObject_name}.h5ad")

