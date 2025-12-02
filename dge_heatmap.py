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

    if "rank_genes_groups" not in adata.uns:
        raise RuntimeError("rank_genes_groups missing. Run dge.py first.")

    r = adata.uns["rank_genes_groups"]
    scanpy_groups = r["names"].dtype.names  # original scanpy groups

    # ------------------------
    # Build dataframe of DE results
    # ------------------------
    rows = []
    for g in scanpy_groups:
        for i in range(len(r["names"][g])):
            rows.append({
                "group": g,
                "gene": r["names"][g][i],
                "score": r["scores"][g][i],
            })

    df = pd.DataFrame(rows)

    # Top N per group
    top_df = (
        df.sort_values("score", ascending=False)
          .groupby("group")
          .head(topn)
    )

    genes = top_df["gene"].unique().tolist()

    # ------------------------
    # Extract expression
    # ------------------------
    if adata.raw is not None:
        X = adata.raw[:, genes].X
    else:
        X = adata[:, genes].X

    if hasattr(X, "toarray"):
        X = X.toarray()

    # z-score across cells
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-9)

    df_expr = pd.DataFrame(X, columns=genes)
    df_expr[groupby] = adata.obs[groupby].values

    # ------------------------
    # Custom sorting for celltype
    # ------------------------
    if groupby == "celltype":
        custom_order = (
            'MG', 'MGPC', 'PR precursors', 'Rod', 'Cones',
            'BC', 'AC', 'HC', 'RGC', 'Microglia_ImmuneCells','RPE','Endothelial', 'Pericytes','Oligocytes','Melanocytes'
        )
        # keep only those present in data
        present_order = [x for x in custom_order if x in df_expr[groupby].unique()]
        df_expr[groupby] = pd.Categorical(df_expr[groupby], categories=present_order, ordered=True)
        df_expr = df_expr.sort_values(groupby)
        ordered_groups = present_order
    else:
        df_expr = df_expr.sort_values(groupby)
        ordered_groups = df_expr[groupby].unique().tolist()

    # Extract heatmap matrix
    heatmap_data = df_expr.drop(columns=[groupby]).to_numpy()

    # ------------------------
    # Plot heatmap
    # ------------------------
    plt.figure(figsize=(12, 16))
    plt.imshow(heatmap_data.T, aspect='auto', cmap='RdBu_r', interpolation='nearest')
    plt.colorbar(label='z-score')

    # gene names along y-axis
    plt.yticks(range(len(genes)), genes, fontsize=4)

    # Correct x-axis labeling: one tick per group
    group_cell_counts = df_expr[groupby].value_counts()[ordered_groups].tolist()
    group_tick_positions = np.cumsum([0] + group_cell_counts[:-1])

    plt.xticks(
        ticks=group_tick_positions,
        labels=ordered_groups,
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

