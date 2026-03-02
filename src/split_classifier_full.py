#!/usr/bin/env python3
# split_classifier_full_fixed2.py
# Fully robust version – no view warnings, handles sparse matrix copy errors,
# and prints diagnostic info to debug empty violin plots.

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import OneClassSVM
from sklearn.metrics import confusion_matrix
from scipy.sparse import issparse, csr_matrix
import warnings
import traceback

warnings.simplefilter("ignore", RuntimeWarning)

def to_dense(X):
    return X.A if issparse(X) else X

def valid_cells(X):
    X_dense = to_dense(X)
    return np.all(np.isfinite(X_dense), axis=1)

def safe_subset_copy(adata, idx):
    """
    Attempt to create a proper copy of adata[idx].
    If the standard .copy() fails (due to sparse matrix bugs),
    manually reconstruct the AnnData.
    """
    try:
        return adata[idx].copy()
    except Exception as e:
        print(f"Standard copy failed: {e}. Attempting manual reconstruction...")
        # Manual reconstruction
        X = adata.X[idx] if issparse(adata.X) else adata.X[idx]
        if issparse(X):
            X = X.copy()  # ensure own data
        else:
            X = np.array(X)
        obs = adata.obs.iloc[idx].copy()
        var = adata.var.copy()  # var is same for all cells
        obsm = {k: v[idx] for k, v in adata.obsm.items()} if hasattr(adata, 'obsm') else {}
        uns = adata.uns  # uns is shared, but we can keep as is (read-only)
        new_adata = sc.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
        return new_adata

def run_control_similarity(sub, condition, hvg_n, max_control_cells=3000, random_state=0):
    """
    sub: AnnData object (copy) containing only one cell type.
    condition: "LD" or "NMDA"
    Returns (ctrl, test) where each is an AnnData with fidelity scores in .obs.
    """
    # Extract control and test cells – make copies to avoid views
    ctrl_all = sub[sub.obs["renamed_samples"] == "Control"].copy()
    test = sub[sub.obs["renamed_samples"] == condition].copy()

    # Remove Inf/NaN Control cells
    valid_mask = valid_cells(ctrl_all.X)
    ctrl_all = ctrl_all[valid_mask].copy()

    if ctrl_all.n_obs < 10 or test.n_obs == 0:
        print(f"   Skipping: not enough valid cells for {condition} (Control {ctrl_all.n_obs}, Test {test.n_obs})")
        return None, None

    # Downsample Control
    if ctrl_all.n_obs > max_control_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(ctrl_all.n_obs, max_control_cells, replace=False)
        ctrl = ctrl_all[idx].copy()
    else:
        ctrl = ctrl_all.copy()

    # Remove zero-variance genes
    keep = np.var(to_dense(ctrl.X), axis=0) > 0
    ctrl = ctrl[:, keep].copy()
    test = test[:, keep].copy()

    # HVG
    hvg_n = min(hvg_n, ctrl.n_vars)
    try:
        sc.pp.highly_variable_genes(ctrl, n_top_genes=hvg_n, flavor="seurat", subset=True)
    except Exception:
        var_idx = np.argsort(np.var(to_dense(ctrl.X), axis=0))[-hvg_n:]
        ctrl = ctrl[:, var_idx].copy()
    # Align test genes
    test = test[:, ctrl.var_names].copy()   # <-- crucial: copy after gene selection

    X_train = to_dense(ctrl.X)
    X_test = to_dense(test.X)

    # Remove NaN/Inf
    X_train = np.nan_to_num(X_train)
    X_test = np.nan_to_num(X_test)

    model = Pipeline([
        ("scaler", StandardScaler()),
        ("ocsvm", OneClassSVM(kernel="rbf", nu=0.05, gamma="scale"))
    ])
    model.fit(X_train)

    f_train = model.decision_function(X_train)
    f_test = model.decision_function(X_test)

    all_f = np.concatenate([f_train, f_test])
    f_scaled = (all_f - all_f.min()) / (all_f.max() - all_f.min() + 1e-9)

    # Now both ctrl and test are proper copies, so we can assign to .obs safely
    ctrl.obs[f"fidelity_{condition}"] = f_scaled[:len(f_train)]
    test.obs[f"fidelity_{condition}"] = f_scaled[len(f_train):]

    y_true = np.array([1]*len(f_train) + [0]*len(f_test))
    y_pred = (f_scaled >= 0.5).astype(int)
    print(f"\nConfusion matrix: Control vs {condition} (subset)")
    print(confusion_matrix(y_true, y_pred))

    return ctrl, test

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--max_control_cells", type=int, default=3000)
    parser.add_argument("--celltype_key", default="celltype",
                        help="Column to split subsets by cell type (e.g., Rod, BC, HC)")
    args = parser.parse_args()

    suffix = f"_{args.suffix}" if args.suffix else ""
    adata = sc.read_h5ad(args.h5ad)

    # Initialize fidelity columns with NaN
    adata.obs["fidelity_LD"] = np.nan
    adata.obs["fidelity_NMDA"] = np.nan

    if args.celltype_key not in adata.obs.columns:
        raise ValueError(f"Cell type key '{args.celltype_key}' not found in adata.obs")

    celltypes = adata.obs[args.celltype_key].unique()
    for ct in celltypes:
        print(f"\nProcessing cell type {ct} ...")

        mask = (adata.obs[args.celltype_key] == ct).values
        if not np.any(mask):
            print(f"Cell type {ct} has no cells, skipping.")
            continue
        indices = np.where(mask)[0]

        # Use safe copy to avoid sparse matrix errors
        try:
            sub = safe_subset_copy(adata, indices)
        except Exception as e:
            print(f"Error creating subset for cell type {ct}: {e}")
            traceback.print_exc()
            print("Skipping this cell type.")
            continue

        # Check that sub contains both Control and at least one condition
        conditions_present = sub.obs["renamed_samples"].unique()
        print(f"   Conditions present: {conditions_present}")
        print(f"   Control count: {(sub.obs['renamed_samples']=='Control').sum()}")
        print(f"   LD count: {(sub.obs['renamed_samples']=='LD').sum()}")
        print(f"   NMDA count: {(sub.obs['renamed_samples']=='NMDA').sum()}")

        try:
            ctrl_ld, ld = run_control_similarity(sub, "LD", args.hvg_genes,
                                                 max_control_cells=args.max_control_cells)
            ctrl_nmda, nmda = run_control_similarity(sub, "NMDA", args.hvg_genes,
                                                     max_control_cells=args.max_control_cells)
        except Exception as e:
            print(f"Error processing cell type {ct}: {e}")
            traceback.print_exc()
            print("Skipping this cell type.")
            continue

        # Assign back to original adata using .loc with .values
        if ctrl_ld is not None:
            adata.obs.loc[ctrl_ld.obs_names, "fidelity_LD"] = ctrl_ld.obs["fidelity_LD"].values
        if ld is not None:
            adata.obs.loc[ld.obs_names, "fidelity_LD"] = ld.obs["fidelity_LD"].values
        if ctrl_nmda is not None:
            adata.obs.loc[ctrl_nmda.obs_names, "fidelity_NMDA"] = ctrl_nmda.obs["fidelity_NMDA"].values
        if nmda is not None:
            adata.obs.loc[nmda.obs_names, "fidelity_NMDA"] = nmda.obs["fidelity_NMDA"].values

    # Diagnostic: count non-NaN fidelity values per condition
    print("\n=== FINAL DIAGNOSTICS ===")
    for cond in ["LD", "NMDA"]:
        col = f"fidelity_{cond}"
        for samp in ["Control", "LD", "NMDA"]:
            mask = adata.obs["renamed_samples"] == samp
            n_nonnan = np.sum(~np.isnan(adata.obs.loc[mask, col]))
            print(f"{col} for {samp}: {n_nonnan} non-NaN values out of {np.sum(mask)} cells")

    # UMAP (if not already present)
    if "X_umap" not in adata.obsm:
        print("Computing UMAP...")
        sc.pp.neighbors(adata, n_jobs=4)
        sc.tl.umap(adata, n_jobs=4)

    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.legend(loc="best", frameon=False)
    plt.savefig(f"umap_all_conditions{suffix}.png", dpi=300)
    plt.close()

    # Violin plot – only if there is data
    groups = ["Control", "LD", "NMDA"]
    data = [
        adata.obs.loc[adata.obs["renamed_samples"] == "Control", "fidelity_LD"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "LD", "fidelity_LD"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "NMDA", "fidelity_NMDA"].dropna().values
    ]
    for i, d in enumerate(data):
        if len(d) == 0:
            print(f"Warning: No data for group {groups[i]}. Violin will be empty for this group.")

    if any(len(d) > 0 for d in data):
        fig, ax = plt.subplots(figsize=(6, 4))
        vp = ax.violinplot(data, positions=range(len(groups)), widths=0.6,
                           showmeans=False, showextrema=False, showmedians=True)
        for b in vp["bodies"]:
            b.set_alpha(0.7)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups)
        ax.set_ylabel("Control similarity (fidelity)")
        ax.set_ylim(0, 1)
        plt.tight_layout()
        plt.savefig(f"violin_fidelity_all_groups{suffix}.png", dpi=300)
        plt.close()
    else:
        print("No fidelity data at all – violin plot not generated.")

    # Save final AnnData
    adata.write(args.output)
    print(f"Saved results to {args.output}")

if __name__ == "__main__":
    main()
