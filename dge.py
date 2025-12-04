#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
    parser.add_argument("-c", "--csv", required=True, help="Output CSV file for full DE table")
    parser.add_argument("-g", "--groupby", required=True, help="Column in adata.obs to group by")
    parser.add_argument("-o", "--out", required=True, help="Output H5AD file with DE results")
    parser.add_argument("-n", "--topn", type=int, default=None, help="Number of top genes per cluster")
    parser.add_argument("-m", "--metric", type=str, default="wilcoxon_score", help="Metric to rank genes by (wilcoxon_score, logfoldchange, pval, pval_adj)")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs:
        raise KeyError(f"Groupby key '{args.groupby}' not found in adata.obs")

    # Fix for log1p metadata
    adata_copy = adata.copy()
    adata_copy.uns["log1p"] = {"base": None}

    # Run DE
    sc.tl.rank_genes_groups(
        adata_copy,
        groupby=args.groupby,
        method="wilcoxon",
        n_genes=None
    )

    # Convert results → dataframe
    r = adata_copy.uns["rank_genes_groups"]
    groups = r["names"].dtype.names

    rows = []
    for g in groups:
        for i in range(len(r["names"][g])):
            rows.append({
                "group": g,
                "gene": r["names"][g][i],
                "wilcoxon_score": r["scores"][g][i],
                "pval": r["pvals"][g][i],
                "pval_adj": r["pvals_adj"][g][i],
                "logfoldchange": r["logfoldchanges"][g][i]
            })

    df = pd.DataFrame(rows)
    df.to_csv(args.csv, index=False)

    # Save H5AD with DE results included
    adata_copy.uns["dge_table"] = df
    adata_copy.write_h5ad(args.out)

    print(f"Saved full DGE table → {args.csv}")
    print(f"Saved DE-containing object → {args.out}")

    # Optional: Save top n genes per cluster based on metric
    if args.topn is not None:
        if args.metric not in df.columns:
            raise ValueError(f"Metric '{args.metric}' not found in DE table columns: {list(df.columns)}")

        top_rows = []
        for g, group_df in df.groupby("group"):
            # Sort by absolute value of the chosen metric, descending
            top_group = group_df.reindex(
                group_df[args.metric].abs().sort_values(ascending=False).index
            ).head(args.topn)
            top_rows.append(top_group)

        df_top = pd.concat(top_rows)

        # Construct new CSV name: original + _{n}{metric}.csv
        base, ext = os.path.splitext(args.csv)
        top_csv = f"{base}_{args.topn}{args.metric}{ext}"
        df_top.to_csv(top_csv, index=False)
        print(f"Saved top {args.topn} genes per cluster by absolute '{args.metric}' → {top_csv}")

if __name__ == "__main__":
    main()

