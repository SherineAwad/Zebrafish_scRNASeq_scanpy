#!/usr/bin/env python3

import pandas as pd
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="DGE: group vs rest (Scanpy results) + heatmap")
    parser.add_argument("-i", "--ih5ad", required=True, help="Input .h5ad file")
    parser.add_argument("-c", "--csv", required=True, help="Output CSV file")
    parser.add_argument("-o", "--oh5ad", default=None, help="Optional output .h5ad")
    parser.add_argument("-g", "--groupby", required=True, help="Groupby column for DGE")
    parser.add_argument("-n", "--topn", type=int, default=20, help="Number of top genes per group to plot in heatmap")
    parser.add_argument("--heatmap", default="heatmap.png", help="Output heatmap PNG filename")

    args = parser.parse_args()

    groupby = args.groupby
    topn = args.topn

    # Load AnnData
    adata = sc.read_h5ad(args.ih5ad)

    # Check for groupby annotation
    if groupby not in adata.obs:
        raise ValueError(f"adata.obs[{groupby}] not found")

    # Make a copy to avoid KeyError from pre-normalized X
    adata_for_dge = adata.copy()
    adata_for_dge.uns['log1p'] = {'base': None}

    # Run Wilcoxon DGE: each group vs all others, using raw data if available
    sc.tl.rank_genes_groups(
        adata_for_dge,
        groupby=groupby,
        method="wilcoxon",
        reference="rest",
        use_raw=True
    )

    # Collect results into a DataFrame exactly as Scanpy outputs
    groups = adata_for_dge.obs[groupby].unique()
    all_results = []

    result = adata_for_dge.uns['rank_genes_groups']
    for group in groups:
        names = result['names'][group]
        scores = result['scores'][group]
        pvals_adj = result['pvals_adj'][group]
        logfoldchanges = result['logfoldchanges'][group]

        df = pd.DataFrame({
            'gene': names,
            'cell_type': group,
            'wilcoxon_score': scores,
            'p_val_adj': pvals_adj,
            'logfoldchange': logfoldchanges,
            'method': 'wilcoxon'
        })
        all_results.append(df)

    df_all = pd.concat(all_results, ignore_index=True)
    df_all.to_csv(args.csv, index=False)
    print(f"Saved DGE results → {args.csv}")

    # ---------------------------------------------------
    # Heatmap of top N genes per group
    # ---------------------------------------------------
    # Collect top genes across all groups
    top_genes = df_all.groupby("cell_type").apply(lambda x: x.nlargest(topn, "wilcoxon_score")).reset_index(drop=True)
    genes_to_plot = top_genes['gene'].unique()

    # Extract expression matrix
    if adata.raw is not None:
        expr = adata.raw[:, genes_to_plot].X
    else:
        expr = adata[:, genes_to_plot].X

    # Convert sparse to dense if needed
    if hasattr(expr, "toarray"):
        expr = expr.toarray()

    # z-score normalization per gene for heatmap
    expr = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-9)

    # Create DataFrame for plotting
    df_expr = pd.DataFrame(expr, index=adata.obs_names, columns=genes_to_plot)
    df_expr[groupby] = adata.obs[groupby]

    # Sort by group
    df_expr = df_expr.sort_values(by=groupby)
    heatmap_data = df_expr.drop(columns=[groupby])

    # Plot heatmap
    plt.figure(figsize=(12, 8))
    plt.imshow(heatmap_data.T, aspect='auto', cmap='RdBu_r', interpolation='nearest')
    plt.colorbar(label='z-score')
    plt.yticks(range(len(genes_to_plot)), genes_to_plot,fontsize=5)
    plt.xticks([])
    plt.title(f"Top {topn} genes per {groupby}")
    plt.tight_layout()
    plt.savefig(args.heatmap, dpi=300)
    plt.close()
    print(f"Saved heatmap → {args.heatmap}")

    # Save AnnData if requested
    if args.oh5ad:
        adata.write(args.oh5ad)
        print(f"Saved AnnData → {args.oh5ad}")

if __name__ == "__main__":
    main()

