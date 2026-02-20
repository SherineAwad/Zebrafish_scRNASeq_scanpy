#!/usr/bin/env python3
# classifier_control_vs_ld_nmda_continuous.py
# Full fix: Control-only training per comparison with graded fidelity

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import OneClassSVM
from sklearn.metrics import confusion_matrix

def run_control_similarity(adata, condition, hvg_n, suffix):
    """
    Train a Control-only One-Class SVM and compute graded fidelity scores
    for Control vs <condition>
    """

    # Subset cells
    train_adata = adata[adata.obs["renamed_samples"] == "Control"].copy()
    test_adata  = adata[adata.obs["renamed_samples"] == condition].copy()

    if train_adata.n_obs < 10:
        raise RuntimeError("Not enough Control cells to train.")

    # -------------------------
    # Remove zero-variance genes (based on Control)
    # -------------------------
    Xc = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
    gene_var = np.var(Xc, axis=0)
    keep = gene_var > 0
    train_adata = train_adata[:, keep].copy()
    test_adata  = test_adata[:, keep].copy()

    # -------------------------
    # HVG selection (Control only)
    # -------------------------
    hvg_n = min(hvg_n, train_adata.n_vars)
    try:
        sc.pp.highly_variable_genes(train_adata, n_top_genes=hvg_n, flavor="seurat", subset=True)
        X_train = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
        X_test  = test_adata.X[:, train_adata.var.highly_variable.values]
    except Exception:
        X_train = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
        X_test  = test_adata.X.toarray() if not isinstance(test_adata.X, np.ndarray) else test_adata.X

    # -------------------------
    # NaN / Inf safety
    # -------------------------
    X_train = np.nan_to_num(X_train, nan=0.0, posinf=1e6, neginf=-1e6)
    X_test  = np.nan_to_num(X_test, nan=0.0, posinf=1e6, neginf=-1e6)

    # -------------------------
    # One-Class SVM to measure similarity to Control
    # -------------------------
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("ocsvm", OneClassSVM(kernel="rbf", nu=0.05, gamma="scale"))
    ])

    model.fit(X_train)

    # Decision function: higher = more similar to Control
    f_train = model.decision_function(X_train)
    f_test  = model.decision_function(X_test)

    # Scale to 0–1 for plotting
    all_f = np.concatenate([f_train, f_test])
    f_scaled = (all_f - all_f.min()) / (all_f.max() - all_f.min() + 1e-9)

    train_adata.obs[f"fidelity_{condition}"] = f_scaled[:len(f_train)]
    test_adata.obs[f"fidelity_{condition}"]  = f_scaled[len(f_train):]

    # -------------------------
    # Confusion matrix (optional threshold=0.5)
    # -------------------------
    y_true = np.array([1]*len(f_train) + [0]*len(f_test))
    y_pred = (f_scaled >= 0.5).astype(int)
    cm = confusion_matrix(y_true, y_pred)
    print(f"\nConfusion matrix: Control vs {condition}")
    print(cm)

    # -------------------------
    # Combine for plotting
    # -------------------------
    combined = train_adata.concatenate(test_adata, index_unique=None)

    # -------------------------
    # UMAP
    # -------------------------
    if "X_umap" not in combined.obsm:
        sc.pp.neighbors(combined)
        sc.tl.umap(combined)

    # UMAP colored by fidelity
    sc.pl.umap(combined, color=f"fidelity_{condition}", cmap="viridis",
               vmin=0, vmax=1, show=False)
    plt.savefig(f"umap_fidelity_Control_vs_{condition}{suffix}.png", dpi=300)
    plt.close()

    # Violin plot
    sc.pl.violin(combined, keys=f"fidelity_{condition}",
                 groupby="renamed_samples", stripplot=False, show=False)
    plt.savefig(f"violin_fidelity_Control_vs_{condition}{suffix}.png", dpi=300)
    plt.close()

    return combined

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--suffix", default="")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    suffix = f"_{args.suffix}" if args.suffix else ""

    # -------------------------
    # Run Control similarity for LD and NMDA separately
    # -------------------------
    combined_ld = run_control_similarity(adata, "LD", args.hvg_genes, suffix)
    combined_nmda = run_control_similarity(adata, "NMDA", args.hvg_genes, suffix)

    # -------------------------
    # Save all output
    # -------------------------
    output_adata = combined_ld.concatenate(combined_nmda,
                                           index_unique=None,
                                           batch_key="comparison",
                                           batch_categories=["Control_vs_LD", "Control_vs_NMDA"])
    output_adata.write(args.output)

if __name__ == "__main__":
    main()
