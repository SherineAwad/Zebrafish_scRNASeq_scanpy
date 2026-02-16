#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def main():
    # ------------------------
    # Parse command-line args
    # ------------------------
    parser = argparse.ArgumentParser(
        description="Compute Pearson correlation between two conditions across cell types and plot heatmap"
    )
    parser.add_argument("input_file", type=str, help="Path to the .h5ad file")
    parser.add_argument("condition_column", type=str, help="Name of the column storing conditions (e.g., renamed_samples)")
    parser.add_argument("condition1", type=str, help="First condition name (e.g., Control)")
    parser.add_argument("condition2", type=str, help="Second condition name (e.g., LD)")
    parser.add_argument("--min_mean_expr", type=float, required=True, help="Minimum mean expression for genes to be included")
    parser.add_argument("output_heatmap", type=str, help="Path to save the heatmap (e.g., heatmap.png)")
    args = parser.parse_args()

    # ------------------------
    # Load the AnnData object
    # ------------------------
    adata = sc.read_h5ad(args.input_file)

    # ------------------------
    # Use raw data if available
    # ------------------------
    if adata.raw is not None:
        X = adata.raw.X
        gene_names = adata.raw.var_names
        print("Using raw counts from .raw")
    else:
        X = adata.X
        gene_names = adata.var_names
        print("No .raw found, using adata.X")

    # Convert to DataFrame for easier handling
    if isinstance(X, np.ndarray):
        df = pd.DataFrame(X, index=adata.obs_names, columns=gene_names)
    else:
        # sparse matrix
        df = pd.DataFrame(X.toarray(), index=adata.obs_names, columns=gene_names)

    # ------------------------
    # Define celltypes in exact order
    # ------------------------
    common_celltypes = ['MG', 'MGPC', 'PR precursors', 'Rod', 'Cones', 'BC', 'AC', 'HC', 'RGC']

    print(f"Using celltypes in order: {common_celltypes}")

    # Check which celltypes actually exist in the data
    available_celltypes = set(adata.obs['celltype'].unique())
    existing_celltypes = [ct for ct in common_celltypes if ct in available_celltypes]
    
    if not existing_celltypes:
        raise ValueError("None of the specified celltypes found in the data!")
    
    print(f"Available celltypes in data: {existing_celltypes}")
    common_celltypes = existing_celltypes

    # ------------------------
    # Compute average expression per cell type for each condition
    # ------------------------
    avg_expr_cond1 = {}
    avg_expr_cond2 = {}

    for ct in common_celltypes:
        # Condition 1
        idx_cond1 = (adata.obs['celltype'] == ct) & (adata.obs[args.condition_column] == args.condition1)
        if idx_cond1.sum() > 0:
            avg_expr_cond1[ct] = df.loc[idx_cond1].mean(axis=0)
        else:
            print(f"Warning: No cells found for {ct} in {args.condition1}")
            continue

        # Condition 2
        idx_cond2 = (adata.obs['celltype'] == ct) & (adata.obs[args.condition_column] == args.condition2)
        if idx_cond2.sum() > 0:
            avg_expr_cond2[ct] = df.loc[idx_cond2].mean(axis=0)
        else:
            print(f"Warning: No cells found for {ct} in {args.condition2}")
            continue

    # Remove celltypes that don't have data in both conditions
    valid_celltypes = []
    for ct in common_celltypes:
        if ct in avg_expr_cond1 and ct in avg_expr_cond2:
            valid_celltypes.append(ct)
    
    if not valid_celltypes:
        raise ValueError("No celltypes found with data in both conditions!")
        
    common_celltypes = valid_celltypes
    print(f"Final celltypes with data in both conditions: {common_celltypes}")

    # ------------------------
    # Filter genes based on expression criteria
    # ------------------------
    # Combine all average expression series to check gene expression across all cell types
    all_avg_expr = list(avg_expr_cond1.values()) + list(avg_expr_cond2.values())
    combined_avg_df = pd.concat(all_avg_expr, axis=1)

    # Calculate maximum mean expression across all cell types for each gene
    max_mean_expr_per_gene = combined_avg_df.max(axis=1)

    # Filter genes: mean expression > min_mean_expr in at least one cell type
    genes_to_keep = max_mean_expr_per_gene > args.min_mean_expr
    filtered_genes = gene_names[genes_to_keep]

    total_genes = len(gene_names)
    filtered_genes_count = len(filtered_genes)
    print(f"Total genes: {total_genes}")
    print(f"Genes with mean expression > {args.min_mean_expr} in at least one cell type: {filtered_genes_count} out of {total_genes}")
    print(f"Genes filtered out: {total_genes - filtered_genes_count}")

    if len(filtered_genes) == 0:
        raise ValueError(f"No genes passed the expression filter (min_mean_expr = {args.min_mean_expr}). Try a lower threshold.")

    # ------------------------
    # Filter the average expression data to only include selected genes
    # ------------------------
    for ct in common_celltypes:
        avg_expr_cond1[ct] = avg_expr_cond1[ct].loc[filtered_genes]
        avg_expr_cond2[ct] = avg_expr_cond2[ct].loc[filtered_genes]

    # ------------------------
    # Initialize correlation matrix
    # ------------------------
    correlations = pd.DataFrame(index=common_celltypes, columns=common_celltypes, dtype=float)

    # ------------------------
    # Compute Pearson correlation between all pairs of cell types
    # ------------------------
    for ct1 in common_celltypes:
        for ct2 in common_celltypes:
            # Correlation between condition1-ct1 and condition2-ct2
            corr_value = avg_expr_cond1[ct1].corr(avg_expr_cond2[ct2])
            correlations.loc[ct1, ct2] = corr_value

    # ------------------------
    # Modify output filename to include min_mean_expr
    # ------------------------
    base_name, ext = os.path.splitext(args.output_heatmap)
    output_filename = f"{base_name}_{args.min_mean_expr}{ext}"

    # ------------------------
    # Plot heatmap with fixed scale bar
    # ------------------------
    plt.figure(figsize=(max(8, len(common_celltypes)), max(6, len(common_celltypes))))
    
    # Calculate the actual range of correlation values
    min_corr = correlations.values.min()
    max_corr = 1.0  # Pearson correlation max is always 1
    
    sns.heatmap(correlations, annot=True, cmap="vlag", vmin=min_corr, vmax=max_corr,
                square=True, fmt='.3f', cbar_kws={'label': 'Pearson correlation', 'shrink': 0.8})
    plt.title(f"Pearson correlation: {args.condition1} vs {args.condition2}\n"
              f"(Rows: {args.condition1}, Columns: {args.condition2})\n"
              f"Genes: {filtered_genes_count} out of {total_genes} (mean expr > {args.min_mean_expr})")
    plt.xlabel(f"Cell types in {args.condition2}")
    plt.ylabel(f"Cell types in {args.condition1}")
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    print(f"Heatmap saved to {output_filename}")

    # Print correlation matrix
    print("\nCorrelation matrix:")
    print(correlations)

if __name__ == "__main__":
    main()
