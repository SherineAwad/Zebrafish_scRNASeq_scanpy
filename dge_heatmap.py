#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Plot heatmap from DGE h5ad")
    parser.add_argument("-i", "--ih5ad", required=True, help="Input DGE .h5ad (from dge.py)")
    parser.add_argument("-g", "--groupby", required=True, help="Same groupby used in dge.py")
    parser.add_argument("-n", "--topn", type=int, default=20, help="Top N genes per group")
    parser.add_argument("--heatmap", default="heatmap.png", help="Output PNG filename")
    args = parser.parse_args()

    groupby = args.groupby
    topn = args.topn

    # Load the h5ad containing rank_genes_groups
    adata = sc.read_h5ad(args.ih5ad)

    # Extract DGE results from Scanpy
    if "rank_genes_groups" not in adata.uns:
        raise RuntimeError("rank_genes_groups missing. Run dge.py first.")

    r = adata.uns["rank_genes_groups"]
    groups = r["names"].dtype.names  # all cluster names

    # Build dataframe of Scanpy DE results
    rows = []
    for g in groups:
        for i in range(len(r["names"][g])):
            rows.append({
                "group": g,
                "gene": r["names"][g][i],
                "score": r["scores"][g][i],
            })

    df = pd.DataFrame(rows)

    # Pick topN per group
    top_df = (
        df.sort_values("score", ascending=False)
          .groupby("group")
          .head(topn)
    )

    genes = top_df["gene"].unique().tolist()

    # Expression matrix from raw if exists
    if adata.raw is not None:
        X = adata.raw[:, genes].X
    else:
        X = adata[:, genes].X

    if hasattr(X, "toarray"):
        X = X.toarray()

    # z-score genes
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)

    df_expr = pd.DataFrame(X, columns=genes)
    df_expr[groupby] = adata.obs[groupby].values

    # Sort cells by group
    df_expr = df_expr.sort_values(by=groupby)
    heatmap_data = df_expr.drop(columns=[groupby])
    row_labels = df_expr[groupby].tolist()

    # Plot heatmap
    plt.figure(figsize=(12, 16))
    plt.imshow(heatmap_data.T, aspect='auto', cmap='RdBu_r', interpolation='nearest')
    plt.colorbar(label='z-score')

    # y-axis = genes
    plt.yticks(range(len(genes)), genes, fontsize=4)

    # x-axis = group labels (cluster or celltype)
    plt.xticks(
        ticks=np.linspace(0, heatmap_data.shape[0], num=len(groups), endpoint=False),
        labels=groups,
        rotation=90,
        fontsize=6
    )

    plt.title(f"Top {topn} genes per {groupby}")
    plt.tight_layout()
    plt.savefig(args.heatmap, dpi=300)
    plt.close()

    print(f"Saved heatmap → {args.heatmap}")


if __name__ == "__main__":
    main()

