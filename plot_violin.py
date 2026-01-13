import argparse
import scanpy as sc
import os
import matplotlib.pyplot as plt

# Argument parsing
parser = argparse.ArgumentParser(description="Plot marker genes and UMAPs per sample")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('markers', help='Path to file with marker gene list (one per line)')
args = parser.parse_args()

myObject = args.myObject
markers = args.markers

# Extract base name
parts = os.path.basename(myObject).split("_")
base = os.path.splitext(os.path.basename(myObject))[0]
parts = base.split("_")
base_prefix = "_".join(parts[:3])
plot_suffix ="Harmony.png"
#plot_suffix="noHarmony.png"

# Load marker genes
with open(markers, 'r') as f:
    marker_genes = [g.strip() for g in f if g.strip()]
print("Loaded marker genes:", marker_genes)

# Read AnnData and ensure it's fully in memory
adata = sc.read_h5ad(myObject)

# Print existing celltypes
print("\n=== Existing celltypes in adata ===")
if 'celltype' in adata.obs.columns:
    unique_celltypes = adata.obs['celltype'].unique().tolist()
    print("Celltypes found:")
    for ct in sorted(unique_celltypes):
        print(f"  '{ct}'")
    print(f"\nTotal: {len(unique_celltypes)} unique celltypes")
else:
    print("No 'celltype' column found in adata.obs")
    print("Available columns:", adata.obs.columns.tolist())

print("\n=== adata.obs columns ===")
print(adata.obs.columns.tolist())
print("\n=== adata.obs head ===")
print(adata.obs.head())

print("\n=== adata.var columns ===")
print(adata.var.columns.tolist())
print("\n=== adata.var head ===")
print(adata.var.head())

# If you want to see keys in adata.uns:
print("\n=== adata.uns keys ===")
print(list(adata.uns.keys()))

# Filter marker genes to those present in the data
present = [g for g in marker_genes if g in adata.var_names]
missing = [g for g in marker_genes if g not in adata.var_names]

if missing:
    print(f"Warning: These genes were not found and will be skipped:\n  {missing}")

# Plot dotplot and violin for valid genes
if present:
    sc.pl.dotplot(adata, var_names=present, groupby='celltype', save=f"_{base_prefix}_dotplot{plot_suffix}", show=False)
else:
    print("No valid marker genes found for plotting.")

plt.rcParams['figure.dpi'] = 600

# VIOLIN plots per gene with rotated x-axis labels
# First, let's check what celltypes actually exist in the data
if 'celltype' in adata.obs.columns:
    # Get actual celltypes in the data
    actual_celltypes = adata.obs['celltype'].unique().tolist()
    print(f"\nActual celltypes in data: {actual_celltypes}")
    
    # Create the order list based on what's actually in the data
    # Use only the celltypes from your desired order that actually exist
    desired_order = ['MG', 'MGPC', 'PR precursors', 'Rod', 'Cones', 'BC', 'AC', 'HC', 'RGC', 'Microglia_Immunecells', 'RPE', 'Melanocyte', 'Endothelial', 'Pericytes', 'Oligodendrocyte']
    
    # Filter to only include celltypes that exist in the data
    celltype_order = [ct for ct in desired_order if ct in actual_celltypes]
    
    # Add any remaining celltypes that weren't in the desired order
    remaining_celltypes = [ct for ct in actual_celltypes if ct not in celltype_order]
    celltype_order.extend(sorted(remaining_celltypes))
    
    print(f"Final plotting order: {celltype_order}")
    
    for gene in present:
        sc.pl.violin(adata, keys=gene, groupby='celltype', save=f"_{base_prefix}_{gene}_{plot_suffix}", 
                    show=False, rotation=90, order=celltype_order)
else:
    print("'celltype' column not found, cannot plot violin plots")
