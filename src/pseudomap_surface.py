import argparse
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.spatial.distance import cdist

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Input AnnData file")
    parser.add_argument("--suffix", required=True, help="Suffix for output files")
    parser.add_argument("--nbins", type=int, default=100, help="Number of bins along pseudotime for density smoothing")
    parser.add_argument("--jitter", type=float, default=0.05, help="Vertical jitter for plotting cells")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load AnnData
    adata = sc.read_h5ad(args.h5ad)
    print("✔ Loaded AnnData")

    # ---------------- DiffMap ----------------
    if "X_diffmap" not in adata.obsm:
        if "X_pca" not in adata.obsm or "connectivities" not in adata.obsp:
            raise RuntimeError("PCA or neighbors missing. Cannot compute DiffMap.")
        sc.tl.diffmap(adata, n_comps=10)
        print("✔ DiffMap computed")
    else:
        print("✔ Using existing DiffMap")

    # ---------------- Select root cell automatically ----------------
    # Use DiffMap if available, otherwise PCA
    if "X_diffmap" in adata.obsm:
        embedding = adata.obsm["X_diffmap"]
    else:
        embedding = adata.obsm["X_pca"]
    
    # Find cell closest to origin (0-vector) in embedding space
    distances = cdist(embedding, np.zeros((1, embedding.shape[1])))
    root_cell = np.argmin(distances)
    adata.uns['iroot'] = root_cell
    print(f"✔ Root cell selected automatically: index {root_cell}")

    # ---------------- Pseudotime ----------------
    sc.tl.dpt(adata, n_dcs=10)
    pseudotime = adata.obs['dpt_pseudotime'].values
    pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())  # normalize 0→1

    # ---------------- Y-axis: jittered cell positions ----------------
    y = np.random.uniform(0, 1, size=pseudotime.shape)  # vertical jitter

    # ---------------- 2D Kernel Density ----------------
    xy = np.vstack([pseudotime, y])
    kde = gaussian_kde(xy)
    # Evaluate KDE on grid
    x_grid = np.linspace(0,1,args.nbins)
    y_grid = np.linspace(0,1,args.nbins)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)

    # ---------------- Plot temporal surface ----------------
    plt.figure(figsize=(12,6))
    plt.imshow(Z, origin='lower', aspect='auto', extent=[0,1,0,1], cmap='viridis')
    plt.xlabel("Pseudotime (0 → 1)")
    plt.ylabel("Cell surface (jittered)")
    plt.title("Temporal pseudo-surface map")
    plt.colorbar(label="Cell density")
    plt.tight_layout()
    plt.savefig(f"pseudomap_surface_{args.suffix}.png", dpi=300)
    plt.close()
    print(f"✔ Temporal pseudo-surface map saved: pseudomap_surface_{args.suffix}.png")

    # ---------------- Save updated AnnData ----------------
    adata.write(f"pseudomap_surface_{args.suffix}.h5ad")
    print(f"✔ Updated AnnData saved with DiffMap + pseudotime: pseudomap_surface_{args.suffix}.h5ad")

if __name__ == "__main__":
    main()

