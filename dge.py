#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-c", "--csv", required=True)
    parser.add_argument("-g", "--groupby", required=True)
    parser.add_argument("-o", "--out", required=True)
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

if __name__ == "__main__":
    main()

