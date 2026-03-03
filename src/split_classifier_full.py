#!/usr/bin/env python3
# split_classifier_full_fixed2.py
# Fully robust version — handles sparse matrix copy errors,
# corrected Control fidelity (no forced 1), and fixed violin plot

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import OneClassSVM
from sklearn.metrics import confusion_matrix
from scipy.sparse import issparse
import warnings
import traceback

warnings.simplefilter("ignore", RuntimeWarning)

def to_dense(X):
    return X.A if issparse(X) else X

def valid_cells(X):
    X_dense = to_dense(X)
    return np.all(np.isfinite(X_dense), axis=1)

def safe_subset_copy(adata, idx):
    try:
        return adata[idx].copy()
    except Exception as e:
        print(f"Standard copy failed: {e}. Attempting manual reconstruction...")
        X = adata.X[idx] if issparse(adata.X) else adata.X[idx]
        X = X.copy() if issparse(X) else np.array(X)
        obs = adata.obs.iloc[idx].copy()
        var = adata.var.copy()
        obsm = {k: v[idx] for k, v in adata.obsm.items()} if hasattr(adata, 'obsm') else {}
        uns = adata.uns
        return sc.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)

def run_control_similarity(sub, condition, hvg_n, max_control_cells=3000, random_state=0):
    ctrl_all = sub[sub.obs["renamed_samples"] == "Control"].copy()
    test = sub[sub.obs["renamed_samples"] == condition].copy()

    valid_mask = valid_cells(ctrl_all.X)
    ctrl_all = ctrl_all[valid_mask].copy()

    if ctrl_all.n_obs < 10 or test.n_obs == 0:
        print(f"Skipping {condition}: not enough cells (Control {ctrl_all.n_obs}, Test {test.n_obs})")
        return None, None

    if ctrl_all.n_obs > max_control_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(ctrl_all.n_obs, max_control_cells, replace=False)
        ctrl = ctrl_all[idx].copy()
    else:
        ctrl = ctrl_all.copy()

    keep = np.var(to_dense(ctrl.X), axis=0) > 0
    ctrl = ctrl[:, keep].copy()
    test = test[:, keep].copy()

    hvg_n = min(hvg_n, ctrl.n_vars)
    try:
        sc.pp.highly_variable_genes(ctrl, n_top_genes=hvg_n, flavor="seurat", subset=True)
    except Exception:
        var_idx = np.argsort(np.var(to_dense(ctrl.X), axis=0))[-hvg_n:]
        ctrl = ctrl[:, var_idx].copy()
    test = test[:, ctrl.var_names].copy()

    X_train = np.nan_to_num(to_dense(ctrl.X))
    X_test = np.nan_to_num(to_dense(test.X))

    model = Pipeline([
        ("scaler", StandardScaler()),
        ("ocsvm", OneClassSVM(kernel="rbf", nu=0.01, gamma="scale"))
    ])
    model.fit(X_train)

    f_train = model.decision_function(X_train)
    f_test = model.decision_function(X_test)

    all_f = np.concatenate([f_train, f_test])
    f_scaled = (all_f - all_f.min()) / (all_f.max() - all_f.min() + 1e-9)

    # Assign real scaled SVM scores — no forced 1
    ctrl.obs[f"fidelity_{condition}"] = f_scaled[:len(f_train)]
    test.obs[f"fidelity_{condition}"] = f_scaled[len(f_train):]

    y_true = np.array([1]*len(f_train) + [0]*len(f_test))
    y_pred = (f_scaled >= 0.5).astype(int)
    print(f"\nConfusion matrix: Control vs {condition}")
    print(confusion_matrix(y_true, y_pred))

    return ctrl, test

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--max_control_cells", type=int, default=3000)
    parser.add_argument("--celltype_key", default="celltype")
    args = parser.parse_args()

    suffix = f"_{args.suffix}" if args.suffix else ""
    adata = sc.read_h5ad(args.h5ad)

    adata.obs["fidelity_LD"] = np.nan
    adata.obs["fidelity_NMDA"] = np.nan

    celltypes = adata.obs[args.celltype_key].unique()
    for ct in celltypes:
        mask = (adata.obs[args.celltype_key] == ct).values
        indices = np.where(mask)[0]
        try:
            sub = safe_subset_copy(adata, indices)
        except Exception as e:
            print(f"Error creating subset for {ct}: {e}")
            traceback.print_exc()
            continue

        ctrl_ld, ld = run_control_similarity(sub, "LD", args.hvg_genes, max_control_cells=args.max_control_cells)
        ctrl_nmda, nmda = run_control_similarity(sub, "NMDA", args.hvg_genes, max_control_cells=args.max_control_cells)

        if ctrl_ld is not None:
            adata.obs.loc[ctrl_ld.obs_names, "fidelity_LD"] = ctrl_ld.obs["fidelity_LD"].values
        if ld is not None:
            adata.obs.loc[ld.obs_names, "fidelity_LD"] = ld.obs["fidelity_LD"].values
        if ctrl_nmda is not None:
            adata.obs.loc[ctrl_nmda.obs_names, "fidelity_NMDA"] = ctrl_nmda.obs["fidelity_NMDA"].values
        if nmda is not None:
            adata.obs.loc[nmda.obs_names, "fidelity_NMDA"] = nmda.obs["fidelity_NMDA"].values

    # Violin plot
    groups = ["Control", "LD", "NMDA"]
    data = [
        adata.obs.loc[adata.obs["renamed_samples"] == "Control", "fidelity_LD"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "LD", "fidelity_LD"].clip(0,1).values,
        adata.obs.loc[adata.obs["renamed_samples"] == "NMDA", "fidelity_NMDA"].clip(0,1).values
    ]

    fig, ax = plt.subplots(figsize=(6,4))
    vp = ax.violinplot(data, positions=range(len(groups)), widths=0.6,
                       showmeans=False, showextrema=False, showmedians=True)
    for b in vp["bodies"]:
        b.set_alpha(0.7)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_ylabel("Control similarity (fidelity)")
    ax.set_ylim(0, 1.2)
    ax.set_yticks(np.arange(0,1.3,0.2))
    plt.tight_layout()
    plt.savefig(f"violin_fidelity_all_groups{suffix}.png", dpi=300)
    plt.close()

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.legend(loc="best", frameon=False)
    plt.savefig(f"umap_all_conditions{suffix}.png", dpi=300)
    plt.close()

    adata.write(args.output)
    print(f"Saved results to {args.output}")

if __name__ == "__main__":
    main()
