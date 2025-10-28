#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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
    # Get celltypes present in both conditions
    # ------------------------
    celltypes_cond1 = set(adata.obs[adata.obs[args.condition_column] == args.condition1]['celltype'].unique())
    celltypes_cond2 = set(adata.obs[adata.obs[args.condition_column] == args.condition2]['celltype'].unique())
    common_celltypes = sorted(list(celltypes_cond1 & celltypes_cond2))
    
    print(f"Cell types in {args.condition1}: {celltypes_cond1}")
    print(f"Cell types in {args.condition2}: {celltypes_cond2}")
    print(f"Common cell types: {common_celltypes}")
    
    if not common_celltypes:
        raise ValueError("No common cell types found between the two conditions!")

    # ------------------------
    # Initialize correlation matrix
    # ------------------------
    correlations = pd.DataFrame(index=common_celltypes, columns=common_celltypes, dtype=float)

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
            avg_expr_cond1[ct] = pd.Series(0, index=gene_names)
        
        # Condition 2
        idx_cond2 = (adata.obs['celltype'] == ct) & (adata.obs[args.condition_column] == args.condition2)
        if idx_cond2.sum() > 0:
            avg_expr_cond2[ct] = df.loc[idx_cond2].mean(axis=0)
        else:
            avg_expr_cond2[ct] = pd.Series(0, index=gene_names)

    # ------------------------
    # Compute Pearson correlation between all pairs of cell types
    # ------------------------
    for ct1 in common_celltypes:
        for ct2 in common_celltypes:
            # Correlation between condition1-ct1 and condition2-ct2
            corr_value = avg_expr_cond1[ct1].corr(avg_expr_cond2[ct2])
            correlations.loc[ct1, ct2] = corr_value

    # ------------------------
    # Plot heatmap
    # ------------------------
    plt.figure(figsize=(max(8, len(common_celltypes)), max(6, len(common_celltypes))))
    sns.heatmap(correlations, annot=True, cmap="vlag", vmin=-1, vmax=1, 
                square=True, fmt='.3f', cbar_kws={'label': 'Pearson correlation'})
    plt.title(f"Pearson correlation: {args.condition1} vs {args.condition2}\n"
              f"(Rows: {args.condition1}, Columns: {args.condition2})")
    plt.xlabel(f"Cell types in {args.condition2}")
    plt.ylabel(f"Cell types in {args.condition1}")
    plt.tight_layout()
    plt.savefig(args.output_heatmap, dpi=300)
    print(f"Heatmap saved to {args.output_heatmap}")
    
    # Print correlation matrix
    print("\nCorrelation matrix:")
    print(correlations)

if __name__ == "__main__":
    main()
