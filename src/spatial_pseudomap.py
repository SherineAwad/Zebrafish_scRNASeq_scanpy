#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description="Pseudo-spatial surface map using marker genes")
    parser.add_argument("--h5ad", required=True, help="Input h5ad file")
    parser.add_argument("--suffix", required=True, help="Suffix for output files")
    parser.add_argument(
        "--nasal_genes", nargs="+", required=True, help="Genes enriched in nasal region"
    )
    parser.add_argument(
        "--temporal_genes", nargs="+", required=True, help="Genes enriched in temporal region"
    )
    parser.add_argument(
        "--ventral_genes", nargs="+", required=True, help="Genes enriched in ventral region"
    )
    parser.add_argument(
        "--dorsal_genes", nargs="+", required=True, help="Genes enriched in dorsal region"
    )
    return parser.parse_args()

def compute_gradient(adata, pos_genes, neg_genes, axis_name):
    # Take intersection with genes actually in adata
    pos_genes = [g.lower() for g in pos_genes if g.lower() in adata.var_names.str.lower()]
    neg_genes = [g.lower() for g in neg_genes if g.lower() in adata.var_names.str.lower()]
    if not pos_genes and not neg_genes:
        raise ValueError(f"No {axis_name} marker genes found in the dataset")

    # Mean expression per cell
    expr_pos = adata[:, pos_genes].X.mean(axis=1) if pos_genes else np.zeros(adata.n_obs)
    expr_neg = adata[:, neg_genes].X.mean(axis=1) if neg_genes else np.zeros(adata.n_obs)

    # Gradient = pos - neg
    gradient = np.array(expr_pos).flatten() - np.array(expr_neg).flatten()

    # Normalize 0-1
    gradient = (gradient - gradient.min()) / (gradient.max() - gradient.min())
    return gradient

def main():
    args = parse_args()

    # Load data
    adata = sc.read_h5ad(args.h5ad)
    print(f"✔ Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes")

    # Ensure gene names lowercase
    adata.var_names = adata.var_names.str.lower()

    # Compute pseudo X/Y coordinates
    adata.obs["pseudo_X"] = compute_gradient(adata, args.nasal_genes, args.temporal_genes, "X-axis")
    adata.obs["pseudo_Y"] = compute_gradient(adata, args.ventral_genes, args.dorsal_genes, "Y-axis")

    print("✔ Computed pseudo-spatial coordinates")

    # Plot pseudo-spatial surface with density
    plt.figure(figsize=(6,6))
    plt.hexbin(
        adata.obs["pseudo_X"],
        adata.obs["pseudo_Y"],
        gridsize=100,
        cmap="viridis",
        mincnt=1
    )
    plt.xlabel("Pseudo Nasal → Temporal")
    plt.ylabel("Pseudo Ventral → Dorsal")
    plt.title("Pseudo-Spatial Surface Map")
    plt.colorbar(label="Cell Density")
    plt.tight_layout()
    plt.savefig(f"pseudo_surface_{args.suffix}.png", dpi=300)
    plt.close()

    # Save AnnData with pseudo coordinates
    adata.write(f"pseudo_surface_{args.suffix}.h5ad")
    print("✔ Saved figure and updated AnnData")

if __name__ == "__main__":
    main()

