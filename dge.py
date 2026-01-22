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
    parser.add_argument("-m", "--metric", type=str, default="wilcoxon_score",
                        help="Metric to rank genes by (wilcoxon_score, logfoldchange, pval, pval_adj)")
    parser.add_argument("--sample1", type=str, help="First sample/group for comparison")
    parser.add_argument("--sample2", type=str, help="Second sample/group for comparison")

    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs:
        raise KeyError(f"Groupby key '{args.groupby}' not found in adata.obs")

    if "renamed_samples" not in adata.obs:
        raise KeyError("Column 'renamed_samples' not found in adata.obs")

    adata_copy = adata.copy()
    adata_copy.uns["log1p"] = {"base": None}

    rows = []

    # ================= ORIGINAL ALL-DGE =================
    sc.tl.rank_genes_groups(
        adata_copy,
        groupby=args.groupby,
        method="wilcoxon",
        n_genes=None
    )

    r = adata_copy.uns["rank_genes_groups"]
    groups = r["names"].dtype.names

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

    # ================= TOP-N OVERWRITE =================
    if args.topn is not None:
        if args.metric not in df.columns:
            raise ValueError(f"Metric '{args.metric}' not found in DE table columns: {list(df.columns)}")

        ascending = args.metric in ["pval", "pval_adj"]

        top_rows = []
        for g, group_df in df.groupby("group"):
            if ascending:
                sorted_group = group_df.sort_values(args.metric, ascending=True)
            else:
                sorted_group = group_df.reindex(
                    group_df[args.metric].abs().sort_values(ascending=False).index
                )
            top_rows.append(sorted_group.head(args.topn))

        df = pd.concat(top_rows)

    # ================= SAVE ONLY ONE CSV =================
    df.to_csv(args.csv, index=False)

    adata_copy.uns["dge_table"] = df
    adata_copy.write_h5ad(args.out)

    print(f"Saved DGE table → {args.csv}")
    print(f"Saved DE-containing object → {args.out}")

    # ================= SAMPLE COMPARISON (NO CSV SAVE) =================
    if args.sample1 and args.sample2:

        comp_rows = []

        for cluster in adata_copy.obs[args.groupby].unique():
            sub = adata_copy[adata_copy.obs[args.groupby] == cluster].copy()

            c1 = (sub.obs["renamed_samples"] == args.sample1).sum()
            c2 = (sub.obs["renamed_samples"] == args.sample2).sum()

            if c1 < 2 or c2 < 2:
                print(f"Skipping cluster {cluster}: insufficient cells ({args.sample1}={c1}, {args.sample2}={c2})")
                continue

            sc.tl.rank_genes_groups(
                sub,
                groupby="renamed_samples",
                groups=[args.sample1],
                reference=args.sample2,
                method="wilcoxon",
                n_genes=None
            )

            r = sub.uns["rank_genes_groups"]
            g = args.sample1

            for i in range(len(r["names"][g])):
                comp_rows.append({
                    "cluster": cluster,
                    "sample1": args.sample1,
                    "gene": r["names"][g][i],
                    "wilcoxon_score": r["scores"][g][i],
                    "pval": r["pvals"][g][i],
                    "pval_adj": r["pvals_adj"][g][i],
                    "logfoldchange": r["logfoldchanges"][g][i],
                    "sample2": args.sample2
                })

        df_comp = pd.DataFrame(comp_rows)
        adata_copy.uns["dge_comparison"] = df_comp


if __name__ == "__main__":
    main()
