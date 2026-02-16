#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def main():
    parser = argparse.ArgumentParser(description="Plot Wilcoxon scores heatmap from DGE h5ad")
    parser.add_argument("-i", "--ih5ad", required=True, help="Input DGE .h5ad (from dge.py)")
    parser.add_argument("-g", "--groupby", required=True, help="Same groupby used in dge.py")
    parser.add_argument("-n", "--topn", type=int, default=20, help="Top N genes per group")
    parser.add_argument("--heatmap", default="wilcoxon_heatmap.png", help="Output PNG filename")
    args = parser.parse_args()

    groupby = args.groupby
    topn = args.topn

    # Load the h5ad containing rank_genes_groups
    adata = sc.read_h5ad(args.ih5ad)

    # Extract DGE results from Scanpy
    if "rank_genes_groups" not in adata.uns:
        raise RuntimeError("rank_genes_groups missing. Run dge.py first.")

    r = adata.uns["rank_genes_groups"]
    de_groups = r["names"].dtype.names

    # Build dataframe of Wilcoxon scores
    rows = []
    for g in de_groups:
        for i in range(len(r["names"][g])):
            rows.append({
                "group": g,
                "gene": r["names"][g][i],
                "score": r["scores"][g][i],
            })
    df = pd.DataFrame(rows)

    # Pick topN per group by score
    top_df = (
        df.sort_values("score", ascending=False)
        .groupby("group")
        .head(topn)
    )

    # Create a pivot table: genes x groups with Wilcoxon scores
    heatmap_df = top_df.pivot_table(
        index='gene',
        columns='group',
        values='score',
        aggfunc='first'
    )

    # Fill NaN with 0 (genes not significant for that group)
    heatmap_df = heatmap_df.fillna(0)

    # Reorder columns by manual order if celltype
    if groupby == "celltype":
        manual_order = [
            "MG", "MGPC", "PR precursors", "Rod", "Cones",
            "BC", "AC", "HC", "RGC",
            "Microglia_ImmuneCells", "Perycites", "Melanocyte",
            "Endothelial", "RPE", "Oligodenrocyte"
        ]
        # Filter to existing columns
        manual_order = [g for g in manual_order if g in heatmap_df.columns]
        heatmap_df = heatmap_df[manual_order]
    
    # Sort leiden clusters numerically
    elif groupby == "leiden":
        # Extract numeric cluster numbers and sort
        cluster_nums = sorted([int(x) for x in heatmap_df.columns if x.isdigit()])
        leiden_order = [str(x) for x in cluster_nums]
        heatmap_df = heatmap_df[leiden_order]

    # Sort genes by max score across groups
    heatmap_df = heatmap_df.loc[heatmap_df.max(axis=1).sort_values(ascending=False).index]

    # Create the EXACT scale from your R code (dark blue to light blue, white, light red to dark red)
    colors = [
        '#2166AC',  # Dark blue
        '#67A9CF',  # Medium blue  
        '#D1E5F0',  # Light blue
        '#FFFFFF',  # White
        '#FDDBC7',  # Light red
        '#EF8A62',  # Medium red
        '#B2182B'   # Dark red
    ]

    custom_cmap = LinearSegmentedColormap.from_list('custom_scale', colors, N=256)

    # Plot heatmap
    plt.figure(figsize=(12, 16))
    plt.imshow(heatmap_df.values, aspect='auto', cmap=custom_cmap, interpolation='nearest')
    plt.colorbar(label='Wilcoxon score')

    # y-axis = genes
    plt.yticks(range(len(heatmap_df)), heatmap_df.index, fontsize=4)

    # x-axis = groups
    plt.xticks(range(len(heatmap_df.columns)), heatmap_df.columns, rotation=90, fontsize=6)

    plt.title(f"Top {topn} genes per {groupby} - Wilcoxon scores")
    plt.tight_layout()
    plt.savefig(args.heatmap, dpi=300)
    plt.close()

    print(f"Saved Wilcoxon scores heatmap → {args.heatmap}")

if __name__ == "__main__":
    main()
