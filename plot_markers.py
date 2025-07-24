import argparse
import scanpy as sc
import os

# Argument parsing
parser = argparse.ArgumentParser(description="Plot marker genes and UMAPs per sample")
parser.add_argument('myObject', help='Path to input AnnData (.h5ad) file')
parser.add_argument('markers', help='Path to file with marker gene list (one per line)')
args = parser.parse_args()

myObject = args.myObject
markers = args.markers

# Extract base name
parts = os.path.basename(myObject).split("_")
base_prefix = parts[0]
plot_suffix = ".png"

# Load marker genes
with open(markers, 'r') as f:
    marker_genes = [g.strip() for g in f if g.strip()]
print("Loaded marker genes:", marker_genes)

# Read AnnData and ensure it's fully in memory
adata = sc.read_h5ad(myObject)
print("=== adata.obs columns ===")
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
    sc.pl.violin(adata, keys=present, groupby='celltype', rotation=90, stripplot=True,
                 save=f"_{base_prefix}_violin{plot_suffix}", show=False)
else:
    print("No valid marker genes found for plotting.")

# UMAP feature plots per gene
for gene in present:
    sc.pl.umap(adata, color=gene, save=f"_{base_prefix}_{gene}_featureplot{plot_suffix}", show=False)

renamed_samples = adata.obs['renamed_samples'].cat.categories

for sample in renamed_samples:
    sample_adata = adata[adata.obs['renamed_samples'] == sample].copy()
    # This removes all categories except the ones actually present in the subset
    sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].cat.remove_unused_categories()
    
    sc.pl.umap(
        sample_adata,
        color='renamed_samples',  # Now only one category remains here
        title=f"Sample: {sample}",
        size=20,
        save=f"_{base_prefix}_{sample}_umap.png",
        show=False
    )

# Plot all samples together (uses default color map)
sc.pl.umap(
    adata,
    color='renamed_samples',
    size=2,
    save=f"_{base_prefix}_all_conditions_umap{plot_suffix}",
    show=False
)







