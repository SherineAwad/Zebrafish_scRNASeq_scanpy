#!/usr/bin/env python
import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="DotPlot with gene group annotations")
    parser.add_argument('-i', '--input', required=True, help='Input .h5ad file')
    args = parser.parse_args()

    # Load AnnData
    adata = sc.read_h5ad(args.input)

    # Determine cluster column
    cluster_col = 'combined_leiden' if 'combined_leiden' in adata.obs else 'leiden'

    # ----------------------------
    # Marker genes by functional group
    # ----------------------------
    markers = {
        'GABAergic': ['gad2', 'gad1b', 'slc6a1b'],
        'Cholinergic': ['chata', 'slc18a3a'],
        'glycinergic': ['slc6a9'],
        'dopaminergic': ['th'],
        #'neuropeptide': ['tac1', 'npy', 'vip']
    }

    # Check which genes exist in the dataset
    all_genes = [g for genes in markers.values() for g in genes]
    present_genes = [g for g in all_genes if g in adata.var_names]
    missing_genes = set(all_genes) - set(present_genes)
    if missing_genes:
        print(f"Warning: These genes are missing in the dataset and will be skipped: {missing_genes}")

    # Create the proper structure for Scanpy dotplot with annotations
    # This format is: {group_name: [list_of_genes]}
    markers_to_plot = {}
    for group, gene_list in markers.items():
        # Filter to only include genes that exist
        existing_genes = [g for g in gene_list if g in adata.var_names]
        if existing_genes:  # Only add group if there are genes present
            markers_to_plot[group] = existing_genes

    if not markers_to_plot:
        print("Error: No marker genes found in the dataset!")
        return

    print(f"Plotting genes with groups: {list(markers_to_plot.keys())}")

    # ----------------------------
    # Sort clusters by expression
    # ----------------------------
    print("Sorting clusters by expression...")
    
    # Get all genes from markers
    all_marker_genes = [gene for genes in markers_to_plot.values() for gene in genes]
    
    # Calculate mean expression per cluster for marker genes
    cluster_means = pd.DataFrame(index=adata.obs[cluster_col].cat.categories)
    
    for gene in all_marker_genes:
        if gene in adata.var_names:
            # Get expression data for this gene
            expr = adata[:, gene].X
            if isinstance(expr, np.ndarray):
                expr_array = expr
            else:
                expr_array = expr.toarray()
            
            # Calculate mean expression per cluster
            for cluster in adata.obs[cluster_col].cat.categories:
                mask = adata.obs[cluster_col] == cluster
                if mask.sum() > 0:
                    cluster_means.loc[cluster, gene] = expr_array[mask].mean()
                else:
                    cluster_means.loc[cluster, gene] = 0
    
    # Calculate total expression per cluster (sum across all marker genes)
    cluster_means = cluster_means.astype(float)
    cluster_totals = cluster_means.sum(axis=1)
    
    # Sort clusters by total expression (highest first)
    sorted_clusters = cluster_totals.sort_values(ascending=False).index.tolist()
    
    # Reorder the cluster categories in adata
    adata.obs[cluster_col] = adata.obs[cluster_col].cat.reorder_categories(sorted_clusters)
    
    print(f"Clusters sorted by expression (highest to lowest): {sorted_clusters}")

    # ----------------------------
    # Ensure figures folder exists
    # ----------------------------
    os.makedirs('figures', exist_ok=True)
    base_prefix = os.path.basename(args.input).replace('.h5ad', '')

    # ----------------------------
    # Plot DotPlot with annotations
    # ----------------------------
    print("Generating DotPlot...")

    # Set figure size
    rcParams['figure.figsize'] = (12, 8)

    # Create the dotplot using the dictionary format
    # This is the key: using var_names=markers_to_plot will show annotations
    dot_plot = sc.pl.dotplot(
        adata,
        var_names=markers_to_plot,  # Dictionary with groups as keys
        groupby=cluster_col,
        dendrogram=False,  # Remove dendrogram as requested
        standard_scale='var',
        show=False,
        return_fig=True
    )

    # Save the figure
    output_path = f"figures/{base_prefix}_dotplot_sorted.png"
    dot_plot.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"DotPlot with sorted clusters saved to {output_path}")

    # Also try the Scanpy save mechanism
    try:
        sc.pl.dotplot(
            adata,
            var_names=markers_to_plot,
            groupby=cluster_col,
            dendrogram=False,
            standard_scale='var',
            save=f"_{base_prefix}_annotated_dotplot_sorted.png",
            show=False
        )
        print(f"Also saved via Scanpy to: figures/dotplot_{base_prefix}_annotated_dotplot_sorted.png")
    except Exception as e:
        print(f"Note: Scanpy direct save gave warning: {e}")

if __name__ == "__main__":
    main()
