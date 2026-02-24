#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="Print and plot counts + sample-centered ratios per cluster")
    parser.add_argument("-i", "--input", required=True, help="Input h5ad file")
    parser.add_argument("-c", "--cluster_col", default="leiden", help="Cluster column in adata.obs")
    parser.add_argument("-s", "--sample_col", default="renamed_samples", help="Sample column in adata.obs")
    args = parser.parse_args()

    # Load AnnData
    adata = sc.read_h5ad(args.input)

    # Checks
    if args.cluster_col not in adata.obs:
        raise KeyError(f"Cluster column '{args.cluster_col}' not found in adata.obs")
    if args.sample_col not in adata.obs:
        raise KeyError(f"Sample column '{args.sample_col}' not found in adata.obs")

    # Print total cells
    total_cells = adata.n_obs
    print(f"\nTotal cells: {total_cells}\n")

    # Create output directory
    outdir = "figures"
    os.makedirs(outdir, exist_ok=True)

    # Count cells per cluster per sample
    counts = pd.crosstab(adata.obs[args.cluster_col], adata.obs[args.sample_col])

    # Sample-centered ratios: fraction of sample in each cluster
    ratios = counts.div(counts.sum(axis=0), axis=1)

    # Combine counts and ratios into one column
    combined = counts.astype(str) + " (" + ratios.round(3).astype(str) + ")"
    combined.index.name = "cluster"
    combined.columns.name = "sample"
    combined.to_csv(os.path.join(outdir, "cluster_sample_counts_and_ratios.csv"))

    print("Cell counts and sample-centered ratios per cluster per sample (count (ratio)):\n")
    print(combined)

    # ========== Plot counts ==========
    plt.figure(figsize=(12,6))
    counts.plot(kind="bar", stacked=True, colormap="tab20", width=0.8)
    plt.ylabel("Cell count")
    plt.xlabel("Cluster (leiden)")
    plt.title("Cell counts per cluster per sample")
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    counts_path = os.path.join(outdir, "cluster_sample_counts.png")
    plt.savefig(counts_path, dpi=150)
    plt.close()
    print(f"\nSaved counts plot → {counts_path}")

    # ========== Plot ratios ==========
    plt.figure(figsize=(12,6))
    ratios.plot(kind="bar", stacked=True, colormap="tab20", width=0.8)
    plt.ylabel("Fraction of sample")
    plt.xlabel("Cluster (leiden)")
    plt.title("Sample-centered cell ratios per cluster")
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    ratios_path = os.path.join(outdir, "cluster_sample_ratios.png")
    plt.savefig(ratios_path, dpi=150)
    plt.close()
    print(f"Saved ratios plot → {ratios_path}")

if __name__ == "__main__":
    main()
