#!/usr/bin/env python3
# classifier_control_only_oneclass.py
# TRUE control-only training using One-Class SVM

import argparse
import scanpy as sc
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import OneClassSVM
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--suffix", default="")
    args = parser.parse_args()

    # ----------------------------
    # Load
    # ----------------------------
    adata = sc.read_h5ad(args.h5ad)

    # ----------------------------
    # Subset RGC assumed already done upstream
    # ----------------------------
    train_adata = adata[adata.obs["renamed_samples"] == "Control"].copy()
    test_adata = adata[adata.obs["renamed_samples"].isin(["LD", "NMDA"])].copy()

    if train_adata.n_obs < 10:
        raise RuntimeError("Not enough Control cells to train.")

    # ----------------------------
    # Remove zero-variance genes (based on Control)
    # ----------------------------
    Xc = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
    gene_var = np.var(Xc, axis=0)
    keep = gene_var > 0

    train_adata = train_adata[:, keep].copy()
    test_adata = test_adata[:, keep].copy()

    # ----------------------------
    # HVGs (Control only)
    # ----------------------------
    hvg_n = min(args.hvg_genes, train_adata.n_vars)
    try:
        sc.pp.highly_variable_genes(
            train_adata,
            n_top_genes=hvg_n,
            flavor="seurat",
            subset=False
        )
        hvgs = train_adata.var["highly_variable"].values
        if hvgs.sum() == 0:
            hvgs[:] = True
    except Exception:
        hvgs = np.ones(train_adata.n_vars, dtype=bool)

    # ----------------------------
    # Prepare matrices
    # ----------------------------
    X_train = train_adata.X[:, hvgs]
    X_test = test_adata.X[:, hvgs]

    if not isinstance(X_train, np.ndarray):
        X_train = X_train.toarray()
    if not isinstance(X_test, np.ndarray):
        X_test = X_test.toarray()

    X_train = np.nan_to_num(X_train)
    X_test = np.nan_to_num(X_test)

    # ----------------------------
    # One-class model
    # ----------------------------
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("ocsvm", OneClassSVM(kernel="rbf", nu=0.05, gamma="scale"))
    ])

    model.fit(X_train)

    # ----------------------------
    # Fidelity scores
    # ----------------------------
    f_train = model.decision_function(X_train)
    f_test = model.decision_function(X_test)

    f_all = np.concatenate([f_train, f_test])
    f_all = (f_all - f_all.min()) / (f_all.max() - f_all.min() + 1e-9)

    train_adata.obs["fidelity_score"] = f_all[:len(f_train)]
    test_adata.obs["fidelity_score"] = f_all[len(f_train):]

    # ----------------------------
    # Combine for plotting
    # ----------------------------
    adata_pred = train_adata.concatenate(test_adata, index_unique=None)

    # ----------------------------
    # UMAP
    # ----------------------------
    if "X_umap" not in adata_pred.obsm:
        sc.pp.neighbors(adata_pred)
        sc.tl.umap(adata_pred)

    suffix = f"_{args.suffix}" if args.suffix else ""

    sc.pl.umap(
        adata_pred,
        color="fidelity_score",
        cmap="viridis",
        vmin=0,
        vmax=1,
        show=False
    )
    plt.savefig(f"umap_fidelity{suffix}.png", dpi=300)
    plt.close()

    sc.pl.umap(
        adata_pred,
        color="renamed_samples",
        show=False
    )
    plt.savefig(f"umap_conditions{suffix}.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.violin(
        adata_pred,
        keys="fidelity_score",
        groupby="renamed_samples",
        stripplot=False,
        show=False
    )
    plt.savefig(f"violin_fidelity{suffix}.png", dpi=300)
    plt.close()

    adata_pred.write(args.output)

if __name__ == "__main__":
    main()

