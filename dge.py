#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

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
    parser.add_argument("-k", "--topk", type=int, help="Top K genes for DotPlot")

    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)

    if args.groupby not in adata.obs:
        raise KeyError(f"Groupby key '{args.groupby}' not found in adata.obs")
    if "renamed_samples" not in adata.obs:
        raise KeyError("Column 'renamed_samples' not found in adata.obs")

    adata_copy = adata.copy()
    adata_copy.uns["log1p"] = {"base": None}

    df = None

    # ================ DECIDE: OVERALL SAMPLE COMPARISON OR PER-CLUSTER ================
    if args.sample1 and args.sample2:
        mask = adata_copy.obs["renamed_samples"].isin([args.sample1, args.sample2])
        sub = adata_copy[mask].copy()

        n1 = (sub.obs["renamed_samples"] == args.sample1).sum()
        n2 = (sub.obs["renamed_samples"] == args.sample2).sum()

        if n1 < 2 or n2 < 2:
            print(f"Skipping sample comparison: insufficient cells ({args.sample1}={n1}, {args.sample2}={n2})")
            return

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
        rows = []
        for i in range(len(r["names"][g])):
            rows.append({
                "group": f"{args.sample1}_vs_{args.sample2}",
                "gene": r["names"][g][i],
                "wilcoxon_score": r["scores"][g][i],
                "pval": r["pvals"][g][i],
                "pval_adj": r["pvals_adj"][g][i],
                "logfoldchange": r["logfoldchanges"][g][i]
            })
        df = pd.DataFrame(rows)

        adata_copy.uns["dge_table"] = df
        print(f"Computed OVERALL sample comparison between {args.sample1} and {args.sample2}")

    else:
        rows = []
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
        adata_copy.uns["dge_table"] = df
        print("Computed PER-CLUSTER DGE")

    # ================= TOP-N OVERWRITE =================
    if df is not None and args.topn is not None:
        if args.metric not in df.columns:
            raise ValueError(f"Metric '{args.metric}' not found in DE table columns: {list(df.columns)}")

        finite_df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchange"])

        ascending = args.metric in ["pval", "pval_adj"]
        top_rows = []

        for g, group_df in finite_df.groupby("group"):
            if ascending:
                sorted_group = group_df.sort_values(args.metric, ascending=True)
            else:
                sorted_group = group_df.reindex(
                    group_df[args.metric].abs().sort_values(ascending=False).index
                )
            top_rows.append(sorted_group.head(args.topn))

        df = pd.concat(top_rows)

    # ================= SAVE CSV AND H5AD =================
    if df is not None:
        df.to_csv(args.csv, index=False)
        adata_copy.write_h5ad(args.out)
        print(f"Saved DGE table → {args.csv}")
        print(f"Saved DE-containing object → {args.out}")

    # ================= DOTPLOT (TOP K ONLY FOR PLOTTING) =================
    if args.topk and df is not None:
        fig_dir = "figures"
        os.makedirs(fig_dir, exist_ok=True)

        if args.metric not in df.columns:
            raise ValueError(f"Metric '{args.metric}' not found in DE table columns: {list(df.columns)}")

        finite_df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchange"])

        ascending = args.metric in ["pval", "pval_adj"]

        if args.sample1 and args.sample2:
            if ascending:
                top_genes = finite_df.sort_values(args.metric, ascending=True).head(args.topk)["gene"].tolist()
            else:
                top_genes = finite_df.reindex(
                    finite_df[args.metric].abs().sort_values(ascending=False).index
                ).head(args.topk)["gene"].tolist()
            title = f"{args.sample1}_vs_{args.sample2}"
        else:
            top_genes = []
            for g, group_df in finite_df.groupby("group"):
                if ascending:
                    top_genes.extend(group_df.sort_values(args.metric, ascending=True).head(args.topk)["gene"].tolist())
                else:
                    sorted_group = group_df.reindex(
                        group_df[args.metric].abs().sort_values(ascending=False).index
                    )
                    top_genes.extend(sorted_group.head(args.topk)["gene"].tolist())
            top_genes = list(dict.fromkeys(top_genes))
            title = "AllDGE_topgenes"

        sc.pl.dotplot(
            adata_copy,
            var_names=top_genes,
            groupby=args.groupby,
            standard_scale="var",
            title=title,
            show=False
        )

        csv_prefix = os.path.splitext(os.path.basename(args.csv))[0]
        fig_path = os.path.join(fig_dir, f"{csv_prefix}.png")
        plt.savefig(fig_path, bbox_inches="tight", dpi=150)
        plt.close()
        print(f"Saved DotPlot → {fig_path}")


if __name__ == "__main__":
    main()

