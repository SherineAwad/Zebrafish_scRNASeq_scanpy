#!/usr/bin/env python3
# split_classifier_fixed_full.py
# FULL FIX — ONE Control model for all predictions

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import IsolationForest
from sklearn.metrics import confusion_matrix

def run_control_model(adata, hvg_n):
    """Train ONE model on Control, predict on ALL cells"""
    
    # Get all Control cells for training
    train_adata = adata[adata.obs["renamed_samples"] == "Control"].copy()
    
    if train_adata.n_obs < 10:
        raise RuntimeError("Not enough Control cells to train.")

    # Remove zero-variance genes
    Xc = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
    keep = np.var(Xc, axis=0) > 0
    train_adata = train_adata[:, keep].copy()
    
    # HVG selection on Control only
    hvg_n = min(hvg_n, train_adata.n_vars)
    
    try:
        sc.pp.highly_variable_genes(
            train_adata,
            n_top_genes=hvg_n,
            flavor="seurat",
            subset=True
        )
        X_train = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X
        
    except Exception:
        X_train = train_adata.X.toarray() if not isinstance(train_adata.X, np.ndarray) else train_adata.X

    X_train = np.nan_to_num(X_train)
    
    # Train ONE model on Control only
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("iforest", IsolationForest(contamination=0.01, random_state=42))
    ])
    
    model.fit(X_train)
    
    # Get scores for ALL cells
    all_scores = []
    all_cells = []
    
    for condition in ["Control", "LD", "NMDA"]:
        condition_adata = adata[adata.obs["renamed_samples"] == condition].copy()
        
        # Align genes with training data
        condition_adata = condition_adata[:, train_adata.var_names].copy()
        X_condition = np.nan_to_num(to_dense(condition_adata.X))
        
        # Get decision scores
        scores = model.decision_function(X_condition)
        all_scores.extend(scores)
        all_cells.extend(condition_adata.obs_names)
    
    # Scale all scores together
    all_scores = np.array(all_scores)
    f_scaled = (all_scores - all_scores.min()) / (all_scores.max() - all_scores.min() + 1e-9)
    
    # Create fidelity column
    adata.obs["fidelity"] = np.nan
    for i, cell_name in enumerate(all_cells):
        adata.obs.loc[cell_name, "fidelity"] = f_scaled[i]
    
    # Confusion matrix for Control vs non-Control
    y_true = np.array([1 if adata.obs.loc[cell, "renamed_samples"] == "Control" else 0 for cell in all_cells])
    y_pred = (f_scaled >= 0.5).astype(int)
    
    print("\nConfusion matrix: Control vs all others")
    print(confusion_matrix(y_true, y_pred))
    
    return adata


def to_dense(X):
    from scipy.sparse import issparse
    return X.A if issparse(X) else X


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--hvg_genes", type=int, default=2000)
    parser.add_argument("--suffix", default="")
    args = parser.parse_args()

    suffix = f"_{args.suffix}" if args.suffix else ""

    adata = sc.read_h5ad(args.h5ad)
    
    # Run ONE model on Control, predict everything
    adata = run_control_model(adata, args.hvg_genes)

    # UMAP (independent)
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.legend(loc="best", frameon=False)
    plt.savefig(f"umap_all_conditions{suffix}.png", dpi=300)
    plt.close()

    # Violin plot - ONE model, ALL groups
    groups = ["Control", "LD", "NMDA"]
    x = np.arange(len(groups))

    fig, ax = plt.subplots(figsize=(6,4))

    data = [
        adata.obs.loc[adata.obs["renamed_samples"] == "Control", "fidelity"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "LD", "fidelity"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "NMDA", "fidelity"].dropna().values
    ]

    vp = ax.violinplot(
        data,
        positions=x,
        widths=0.6,
        showmeans=False,
        showextrema=False,
        showmedians=True
    )

    for body in vp["bodies"]:
        body.set_alpha(0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylabel("Control similarity (fidelity)")
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    plt.tight_layout()
    plt.savefig(f"violin_fidelity_all_groups{suffix}.png", dpi=300)
    plt.close()

    adata.write(args.output)
    print(f"Saved results to {args.output}")


if __name__ == "__main__":
    main()
