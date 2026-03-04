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

def run_control_model(adata, hvg_n, contamination, n_estimators, max_samples, max_control_cells, random_state):
    """Train ONE model on Control, predict on ALL cells"""
    
    # Get all Control cells for training
    ctrl_idx = np.where(adata.obs["renamed_samples"] == "Control")[0]
    train_adata = safe_copy_adata(adata, ctrl_idx)
    
    if train_adata.n_obs < 10:
        raise RuntimeError("Not enough Control cells to train.")
    
    # Subsample Control cells if too many
    if train_adata.n_obs > max_control_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(train_adata.n_obs, max_control_cells, replace=False)
        train_adata = train_adata[idx].copy()

    # Remove zero-variance genes
    Xc = to_dense(train_adata.X)
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
    except Exception:
        var_idx = np.argsort(np.var(to_dense(train_adata.X), axis=0))[-hvg_n:]
        train_adata = train_adata[:, var_idx].copy()
    
    X_train = np.nan_to_num(to_dense(train_adata.X))
    
    # Train ONE model on Control only
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("iforest", IsolationForest(
            contamination=contamination,
            n_estimators=n_estimators,
            max_samples=max_samples,
            random_state=random_state
        ))
    ])
    
    model.fit(X_train)
    
    # Get scores for ALL cells
    all_scores = []
    all_cells = []
    
    # Score Control cells (training data) - THESE SHOULD BE HIGH
    train_scores = model.decision_function(X_train)
    all_scores.extend(train_scores)
    all_cells.extend(train_adata.obs_names.tolist())
    
    # Score LD cells
    ld_idx = np.where(adata.obs["renamed_samples"] == "LD")[0]
    if len(ld_idx) > 0:
        ld_adata = safe_copy_adata(adata, ld_idx)
        ld_adata = ld_adata[:, train_adata.var_names].copy()
        X_ld = np.nan_to_num(to_dense(ld_adata.X))
        ld_scores = model.decision_function(X_ld)
        all_scores.extend(ld_scores)
        all_cells.extend(ld_adata.obs_names.tolist())
    
    # Score NMDA cells
    nmda_idx = np.where(adata.obs["renamed_samples"] == "NMDA")[0]
    if len(nmda_idx) > 0:
        nmda_adata = safe_copy_adata(adata, nmda_idx)
        nmda_adata = nmda_adata[:, train_adata.var_names].copy()
        X_nmda = np.nan_to_num(to_dense(nmda_adata.X))
        nmda_scores = model.decision_function(X_nmda)
        all_scores.extend(nmda_scores)
        all_cells.extend(nmda_adata.obs_names.tolist())
    
    # Scale all scores together
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
    parser = argparse.ArgumentParser(description="Isolation Forest model to score cell similarity to Control")
    
    # Input/output parameters
    parser.add_argument("--h5ad", required=True, help="Input AnnData file")
    parser.add_argument("--output", required=True, help="Output AnnData file")
    parser.add_argument("--suffix", default="", help="Suffix for output plots")
    
    # Gene selection parameters
    parser.add_argument("--hvg_genes", type=int, default=3000, 
                       help="Number of highly variable genes to use (default: 3000)")
    
    # Isolation Forest parameters
    parser.add_argument("--contamination", type=float, default=0.01,
                       help="Expected proportion of outliers in Control (default: 0.01)")
    parser.add_argument("--n_estimators", type=int, default=100,
                       help="Number of trees in Isolation Forest (default: 100)")
    parser.add_argument("--max_samples", type=int, default='auto',
                   help="Number of samples per tree (default: 'auto')")
    # Control cell subsampling
    parser.add_argument("--max_control_cells", type=int, default=3000,
                       help="Maximum number of Control cells to use for training (default: 3000)")
    
    # Reproducibility
    parser.add_argument("--random_state", type=int, default=42,
                       help="Random seed for reproducibility (default: 42)")
    
    args = parser.parse_args()

    suffix = f"_{args.suffix}" if args.suffix else ""

    print("Loading data...")
    adata = sc.read_h5ad(args.h5ad)
    
    # Run model with all parameters
    print("Running control model...")
    print(f"Parameters: hvg_genes={args.hvg_genes}, contamination={args.contamination}, "
          f"n_estimators={args.n_estimators}, max_samples={args.max_samples}, "
          f"max_control_cells={args.max_control_cells}, random_state={args.random_state}")
    
    adata = run_control_model(
        adata, 
        hvg_n=args.hvg_genes,
        contamination=args.contamination,
        n_estimators=args.n_estimators,
        max_samples=args.max_samples,
        max_control_cells=args.max_control_cells,
        random_state=args.random_state
    )

    # UMAP
    if "X_umap" not in adata.obsm:
        print("Computing UMAP...")
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    sc.pl.umap(adata, color="renamed_samples", show=False)
    plt.legend(loc="best", frameon=False)
    plt.savefig(f"umap_all_conditions{suffix}.png", dpi=300)
    plt.close()

    # Violin plot - NO median numbers
    print("Creating violin plot...")
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

    # Save
    print(f"Saving results to {args.output}")
    adata.write(args.output)
    print("Done!")

if __name__ == "__main__":
    main()
