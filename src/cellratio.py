#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Print cell counts and ratios per leiden cluster")
    parser.add_argument(
        "-i", "--input", required=True, help="Input h5ad file"
    )
    parser.add_argument(
        "-c", "--cluster_col", default="leiden", help="Column in adata.obs with cluster assignments"
    )
    args = parser.parse_args()

    # Load the AnnData object
    adata = sc.read_h5ad(args.input)

    if args.cluster_col not in adata.obs:
        raise KeyError(f"Cluster column '{args.cluster_col}' not found in adata.obs")

    # Count cells per cluster
    counts = adata.obs[args.cluster_col].value_counts().sort_index()
    total_cells = len(adata)

    # Compute ratios
    ratios = counts / total_cells

    # Combine into a DataFrame for pretty printing
    df = pd.DataFrame({
        "cluster": counts.index,
        "cell_count": counts.values,
        "cell_ratio": ratios.values
    })

    print(f"\nTotal cells: {total_cells}\n")
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
