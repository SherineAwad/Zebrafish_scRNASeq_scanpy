#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--suffix", required=True)
    return p.parse_args()

def main():
    args = parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # -------- checks --------
    if "X_pca" not in adata.obsm:
        raise RuntimeError("PCA missing")
    if "connectivities" not in adata.obsp:
        raise RuntimeError("Neighbors missing")
    if "renamed_samples" not in adata.obs:
        raise RuntimeError("Condition column 'renamed_samples' missing")

    print("✔ Using existing PCA + neighbors")

    # -------- DiffMap --------
    sc.tl.diffmap(adata, n_comps=10)

    # ensure categorical
    adata.obs["renamed_samples"] = adata.obs["renamed_samples"].astype("category")

    # -------- DiffMap colored by condition --------
    sc.pl.diffmap(
        adata,
        components=["1,2"],
        color="renamed_samples",
        show=False
    )
    plt.savefig(
        f"diffmap_{args.suffix}_by_condition.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # -------- DC1 distribution by condition --------
    adata.obs["DC1"] = adata.obsm["X_diffmap"][:, 0]

    sc.pl.violin(
        adata,
        keys="DC1",
        groupby="renamed_samples",
        show=False
    )
    plt.savefig(
        f"diffmap_{args.suffix}_DC1_by_condition.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # -------- save object --------
    adata.write(f"diffmap_{args.suffix}.h5ad")

    print("DONE")
    print("Saved:")
    print(f" diffmap_{args.suffix}_by_condition.png")
    print(f" diffmap_{args.suffix}_DC1_by_condition.png")
    print(f" diffmap_{args.suffix}.h5ad")

if __name__ == "__main__":
    main()

