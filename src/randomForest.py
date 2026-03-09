#!/usr/bin/env python3
# split_classifier_fixed_full.py
# FIXED Random Forest version - NO synthetic data

import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from scipy.sparse import issparse

def to_dense(X):
    return X.toarray() if issparse(X) else X

def safe_copy_adata(adata, idx):
    """Safely copy AnnData with sparse matrix handling"""
    try:
        return adata[idx].copy()
    except Exception:
        X = adata.X[idx]
        if issparse(X):
            X = X.copy()
        else:
            X = np.array(X)
        
        obs = adata.obs.iloc[idx].copy()
        var = adata.var.copy()
        
        obsm = {}
        if hasattr(adata, 'obsm'):
            for key, value in adata.obsm.items():
                try:
                    obsm[key] = value[idx]
                except:
                    obsm[key] = value
        
        return sc.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=adata.uns)

def run_control_model(adata, hvg_n, n_estimators, max_samples, max_control_cells, random_state, min_samples_leaf):
    """Train Random Forest using LD and NMDA as negative class"""
    
    # Get Control cells for training
    ctrl_idx = np.where(adata.obs["renamed_samples"] == "Control")[0]
    train_ctrl = safe_copy_adata(adata, ctrl_idx)
    
    # Get LD cells as negative class
    ld_idx = np.where(adata.obs["renamed_samples"] == "LD")[0]
    train_ld = safe_copy_adata(adata, ld_idx) if len(ld_idx) > 0 else None
    
    if train_ctrl.n_obs < 10 or train_ld is None or train_ld.n_obs < 10:
        raise RuntimeError("Need both Control and LD cells for training.")
    
    # Subsample if too many
    if train_ctrl.n_obs > max_control_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(train_ctrl.n_obs, max_control_cells, replace=False)
        train_ctrl = train_ctrl[idx].copy()
    
    if train_ld.n_obs > max_control_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(train_ld.n_obs, max_control_cells, replace=False)
        train_ld = train_ld[idx].copy()

    # Remove zero-variance genes
    Xc = to_dense(train_ctrl.X)
    keep = np.var(Xc, axis=0) > 0
    train_ctrl = train_ctrl[:, keep].copy()
    train_ld = train_ld[:, keep].copy()
    
    # HVG selection on Control only
    hvg_n = min(hvg_n, train_ctrl.n_vars)
    
    try:
        sc.pp.highly_variable_genes(
            train_ctrl,
            n_top_genes=hvg_n,
            flavor="seurat",
            subset=True
        )
    except Exception:
        var_idx = np.argsort(np.var(to_dense(train_ctrl.X), axis=0))[-hvg_n:]
        train_ctrl = train_ctrl[:, var_idx].copy()
    
    # Align LD with same genes
    train_ld = train_ld[:, train_ctrl.var_names].copy()
    
    X_ctrl = np.nan_to_num(to_dense(train_ctrl.X))
    X_ld = np.nan_to_num(to_dense(train_ld.X))
    
    # Combine Control (1) and LD (0) for training
    X_train = np.vstack([X_ctrl, X_ld])
    y_train = np.hstack([np.ones(len(X_ctrl)), np.zeros(len(X_ld))])
    
    # Train Random Forest classifier
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("rf", RandomForestClassifier(
            n_estimators=n_estimators,
            max_samples=max_samples,
            min_samples_leaf=min_samples_leaf,
            random_state=random_state,
            n_jobs=-1,
            class_weight='balanced'
        ))
    ])
    
    model.fit(X_train, y_train)
    
    # Get scores for ALL cells
    all_scores = []
    all_cells = []
    
    # Score all Control cells
    all_scores.extend(model.predict_proba(X_ctrl)[:, 1])
    all_cells.extend(train_ctrl.obs_names.tolist())
    
    # Score LD cells
    all_scores.extend(model.predict_proba(X_ld)[:, 1])
    all_cells.extend(train_ld.obs_names.tolist())
    
    # Score NMDA cells
    nmda_idx = np.where(adata.obs["renamed_samples"] == "NMDA")[0]
    if len(nmda_idx) > 0:
        nmda_adata = safe_copy_adata(adata, nmda_idx)
        nmda_adata = nmda_adata[:, train_ctrl.var_names].copy()
        X_nmda = np.nan_to_num(to_dense(nmda_adata.X))
        nmda_probs = model.predict_proba(X_nmda)[:, 1]
        all_scores.extend(nmda_probs)
        all_cells.extend(nmda_adata.obs_names.tolist())
    
    # Scale scores
    all_scores = np.array(all_scores)
    f_scaled = (all_scores - all_scores.min()) / (all_scores.max() - all_scores.min() + 1e-9)
    
    # Create fidelity column
    adata.obs["fidelity"] = np.nan
    for i, cell_name in enumerate(all_cells):
        if cell_name in adata.obs_names:
            adata.obs.loc[cell_name, "fidelity"] = f_scaled[i]
    
    # Confusion matrix
    y_true = []
    y_pred = []
    for i, cell_name in enumerate(all_cells):
        if cell_name in adata.obs_names:
            true_label = 1 if adata.obs.loc[cell_name, "renamed_samples"] == "Control" else 0
            y_true.append(true_label)
            y_pred.append(1 if f_scaled[i] >= 0.5 else 0)
    
    print("\nConfusion matrix: Control vs all others")
    print(confusion_matrix(y_true, y_pred))
    
    return adata

def main():
    parser = argparse.ArgumentParser(description="Random Forest model to score cell similarity to Control")
    
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--suffix", default="")
    parser.add_argument("--hvg_genes", type=int, default=3000)
    parser.add_argument("--n_estimators", type=int, default=100)
    parser.add_argument("--max_samples", type=int, default=None)
    parser.add_argument("--min_samples_leaf", type=int, default=1)
    parser.add_argument("--max_control_cells", type=int, default=3000)
    parser.add_argument("--random_state", type=int, default=42)
    
    args = parser.parse_args()

    suffix = f"_{args.suffix}" if args.suffix else ""

    print("Loading data...")
    adata = sc.read_h5ad(args.h5ad)
    
    print("Running Random Forest...")
    print(f"Parameters: hvg_genes={args.hvg_genes}, n_estimators={args.n_estimators}, "
          f"max_samples={args.max_samples}, min_samples_leaf={args.min_samples_leaf}, "
          f"max_control_cells={args.max_control_cells}, random_state={args.random_state}")
    
    adata = run_control_model(
        adata, 
        hvg_n=args.hvg_genes,
        n_estimators=args.n_estimators,
        max_samples=args.max_samples,
        min_samples_leaf=args.min_samples_leaf,
        max_control_cells=args.max_control_cells,
        random_state=args.random_state
    )

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.legend(loc="best", frameon=False)
    plt.savefig(f"umap_all_conditions{suffix}.png", dpi=300)
    plt.close()

    # Violin plot
    groups = ["Control", "LD", "NMDA"]
    x = np.arange(len(groups))

    fig, ax = plt.subplots(figsize=(6,4))

    data = [
        adata.obs.loc[adata.obs["renamed_samples"] == "Control", "fidelity"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "LD", "fidelity"].dropna().values,
        adata.obs.loc[adata.obs["renamed_samples"] == "NMDA", "fidelity"].dropna().values
    ]

    vp = ax.violinplot(data, positions=x, widths=0.6,
                       showmeans=False, showextrema=False, showmedians=True)
    for b in vp["bodies"]:
        b.set_alpha(0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylabel("Control similarity (fidelity)")
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    plt.tight_layout()
    plt.savefig(f"violin_fidelity_all_groups{suffix}.png", dpi=300)
    plt.close()

    adata.write(args.output)
    print("Done!")

if __name__ == "__main__":
    main()
