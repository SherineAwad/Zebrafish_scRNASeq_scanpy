#!/usr/bin/env python3
# split_classifier.py
# FULL FIX — calculations untouched
# UMAP: one plot with Control/LD/NMDA, legend inside
# Violin: one violin per group, shared X/Y axes

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import OneClassSVM
from sklearn.metrics import confusion_matrix

def run_control_similarity(adata, condition, hvg_n):
    train_adata = adata[adata.obs["renamed_samples"] == "Control"].copy()
    test_adata  = adata[adata.obs["renamed_samples"] == condition].copy()

    if train_adata.n_obs < 10:
        raise RuntimeError("Not enough Control cells to train.")

    Xc = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
    keep = np.var(Xc, axis=0) > 0
    train_adata = train_adata[:, keep].copy()
    test_adata  = test_adata[:, keep].copy()

    hvg_n = min(hvg_n, train_adata.n_vars)
    try:
        sc.pp.highly_variable_genes(train_adata, n_top_genes=hvg_n, flavor="seurat", subset=True)
        X_train = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
        X_test  = test_adata.X[:, train_adata.var.highly_variable.values]
    except Exception:
        X_train = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
        X_test  = test_adata.X.toarray() if not isinstance(test_adata.X, np.ndarray) else test_adata.X

    X_train = np.nan_to_num(X_train)
    X_test  = np.nan_to_num(X_test)

    model = Pipeline([("scaler", StandardScaler()),
                      ("ocsvm", OneClassSVM(kernel="rbf", nu=0.05, gamma="scale"))])
    model.fit(X_train)

    f_train = model.decision_function(X_train)
    f_test  = model.decision_function(X_test)

    all_f = np.concatenate([f_train, f_test])
    f_scaled = (all_f - all_f.min()) / (all_f.max() - all_f.min() + 1e-9)

    train_adata.obs[f"fidelity_{condition}"] = f_scaled[:len(f_train)]
    test_adata.obs[f"fidelity_{condition}"]  = f_scaled[len(f_train):]

    y_true = np.array([1]*len(f_train) + [0]*len(f_test))
    y_pred = (f_scaled >= 0.5).astype(int)
    print(f"\nConfusion matrix: Control vs {condition}")
    print(confusion_matrix(y_true, y_pred))

    return train_adata, test_adata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--suffix", default="")
    args = parser.parse_args()

    suffix = f"_{args.suffix}" if args.suffix else ""
    adata = sc.read_h5ad(args.h5ad)

    ctrl_ld, ld = run_control_similarity(adata, "LD", args.hvg_genes)
    ctrl_nmda, nmda = run_control_similarity(adata, "NMDA", args.hvg_genes)

    # Assign fidelities to adata
    adata.obs["fidelity_LD"] = np.nan
    adata.obs["fidelity_NMDA"] = np.nan

    adata.obs.loc[ctrl_ld.obs_names, "fidelity_LD"] = ctrl_ld.obs["fidelity_LD"]
    adata.obs.loc[ld.obs_names, "fidelity_LD"] = ld.obs["fidelity_LD"]

    adata.obs.loc[ctrl_nmda.obs_names, "fidelity_NMDA"] = ctrl_nmda.obs["fidelity_NMDA"]
    adata.obs.loc[nmda.obs_names, "fidelity_NMDA"] = nmda.obs["fidelity_NMDA"]

    # === UMAP: one plot for Control/LD/NMDA, legend inside ===
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.legend(loc="best", frameon=False)  # legend inside plot
    plt.savefig(f"umap_all_conditions{suffix}.png", dpi=300)
    plt.close()

    # === VIOLIN: one violin per group ===
    groups = ["Control", "LD", "NMDA"]
    x = np.arange(len(groups))
    fig, ax = plt.subplots(figsize=(6,4))

    data = [
        adata.obs.loc[adata.obs["renamed_samples"] == "Control", "fidelity_LD"].values,  # Control
        adata.obs.loc[adata.obs["renamed_samples"] == "LD", "fidelity_LD"].values,       # LD
        adata.obs.loc[adata.obs["renamed_samples"] == "NMDA", "fidelity_NMDA"].values    # NMDA
    ]

    vp = ax.violinplot(data, positions=x, widths=0.6, showmeans=False, showextrema=False, showmedians=True)
    for b in vp["bodies"]:
        b.set_alpha(0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylabel("Control similarity (fidelity)")
    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(f"violin_fidelity_all_groups{suffix}.png", dpi=300)
    plt.close()

    adata.write(args.output)

if __name__ == "__main__":
    main()
