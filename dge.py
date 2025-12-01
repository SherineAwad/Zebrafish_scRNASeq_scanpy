#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="DGE: celltype vs rest (Scanpy results only)")
    parser.add_argument("-i", "--ih5ad", required=True, help="Input .h5ad file")
    parser.add_argument("-c", "--csv", required=True, help="Output CSV file")
    parser.add_argument("-o", "--oh5ad", default=None, help="Optional output .h5ad")
    parser.add_argument("-g", "--groupby", required=True, help="groubby" )

    args = parser.parse_args()

    groupby = args.groupby

    # Load AnnData
    adata = sc.read_h5ad(args.ih5ad)

    # Check for groupby annotation
    if groupby not in adata.obs:
        raise ValueError(f"adata.obs[{groupby}] not found")

    # Make a copy to avoid KeyError from pre-normalized X
    adata_for_dge = adata.copy()
    # This line prevents the KeyError in Scanpy
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

    # Save AnnData if requested
    if args.oh5ad:
        adata.write(args.oh5ad)
        print(f"Saved AnnData → {args.oh5ad}")

if __name__ == "__main__":
    main()

