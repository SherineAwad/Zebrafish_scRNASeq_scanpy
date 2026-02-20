#!/usr/bin/env python3
# classifier.py
# FULL FIXED + plots variation for homogeneous subsets
# Arguments / classifier / plotting logic unchanged

import argparse
import scanpy as sc
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    args = parser.parse_args()

    # ----------------------------
    # Load
    # ----------------------------
    adata = sc.read_h5ad(args.h5ad)

    # ----------------------------
    # Keep only Control / LD / NMDA
    # ----------------------------
    adata = adata[adata.obs["renamed_samples"].isin(["Control","LD","NMDA"])].copy()

    # ----------------------------
    # Labels: Control=1, LD/NMDA=0
    # ----------------------------
    y = adata.obs["renamed_samples"].apply(lambda x: 1 if x=="Control" else 0).values
    if len(np.unique(y)) != 2:
        raise RuntimeError("Classifier needs both Control and LD/NMDA cells")

    # ----------------------------
    # Normalize / log1p if needed
    # ----------------------------
    log_needed = True
    if "log1p_total_counts" in adata.obs.columns:
        try:
            if 'log1p' in adata.uns and adata.uns['log1p'].get('base', None) is not None:
                log_needed = False
        except Exception:
            log_needed = True

    if log_needed:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # ----------------------------
    # HVG: remove genes with zero variance
    # ----------------------------
    X_dense = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
    gene_var = np.var(X_dense, axis=0)
    nonzero_var_idx = gene_var > 0

    if nonzero_var_idx.sum() == 0:
        hvgs = np.arange(adata.n_vars)
    else:
        adata = adata[:, nonzero_var_idx].copy()
        hvg_n = min(args.hvg_genes, adata.n_vars)
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg_n, flavor="seurat", subset=False)
        hvgs = adata.var["highly_variable"].values

    # ----------------------------
    # Prepare X and fix NaN/inf
    # ----------------------------
    X = adata.X[:, hvgs]
    if not isinstance(X, np.ndarray):
        X = X.toarray()

    X = np.nan_to_num(X, nan=0.0, posinf=1e6, neginf=-1e6)

    # ----------------------------
    # Classifier
    # ----------------------------
    clf = Pipeline([
        ("scaler", StandardScaler(with_mean=False)),
        ("lr", LogisticRegression(max_iter=2000, solver="lbfgs"))
    ])
    clf.fit(X, y)

    # ----------------------------
    # Fidelity score
    # ----------------------------
    fidelity = clf.predict_proba(X)[:, 1]
    adata.obs["fidelity_score"] = fidelity

    # ----------------------------
    # Add slight jitter for plotting if all identical
    # ----------------------------
    if np.all(fidelity == fidelity[0]):
        jitter = np.random.normal(0, 0.01, size=fidelity.shape)
        adata.obs["fidelity_score_plot"] = np.clip(fidelity + jitter, 0, 1)
    else:
        adata.obs["fidelity_score_plot"] = fidelity

    # ----------------------------
    # UMAP if missing
    # ----------------------------
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    # ----------------------------
    # Plots
    # ----------------------------
    # UMAP colored by fidelity (with jitter for visualization)
    sc.pl.umap(adata, color="fidelity_score_plot", cmap="viridis", vmin=0, vmax=1, show=False)
    plt.savefig("umap_fidelity.png", dpi=300)
    plt.close()

    # UMAP colored by sample
    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.savefig("umap_conditions.png", dpi=300)
    plt.close()

    # Violin plot
    sc.pl.violin(adata, keys="fidelity_score_plot", groupby="renamed_samples",
                 stripplot=False, jitter=False, show=False)
    plt.savefig("violin_fidelity.png", dpi=300)
    plt.close()

    # ----------------------------
    # Save updated h5ad
    # ----------------------------
    adata.write(args.output)

if __name__ == "__main__":
    main()

