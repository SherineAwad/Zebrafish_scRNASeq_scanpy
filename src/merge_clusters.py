#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--annotations', required=True)
    args = parser.parse_args()

    # Read input AnnData
    adata = sc.read_h5ad(args.input)

    # Read annotations file and build mapping
    mapping = {}
    with open(args.annotations, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split(',')
                if len(parts) == 2:
                    old_cluster, new_label = parts[0].strip(), parts[1].strip()
                    mapping[old_cluster] = new_label

    # Apply mapping to leiden column
    leiden_str = adata.obs['leiden'].astype(str).copy()
    adata.obs['combined_leiden'] = leiden_str.map(lambda x: mapping.get(x, x))

    # Save updated AnnData
    adata.write_h5ad(args.output)

    base = args.output.replace('.h5ad', '')

    # -------------------------------
    # Ensure figures/ folder exists
    # -------------------------------
    os.makedirs('figures', exist_ok=True)

    # 1. UMAP colored by sample
    sc.pl.umap(adata, color='renamed_samples', size=80, save=f'{base}_colored_by_sample.png', show=False)

    # 2. UMAP colored by combined_leiden
    sc.pl.umap(adata, color='combined_leiden', size=80, save=f'{base}_combined_leiden.png', show=False)

    # 3. Separate UMAP for each sample with combined_leiden
    for sample in adata.obs['renamed_samples'].unique():
        mask = adata.obs['renamed_samples'] == sample
        adata_sample = adata[mask].copy()
        sc.pl.umap(
            adata_sample,
            color='combined_leiden',
            title=f'Sample: {sample}',
            size=80,
            save=f'{base}_{sample}.png',
            show=False
        )

    # 4. Separate UMAP for each sample WITHOUT leiden, colored distinctly per sample
    plt.rcParams['figure.figsize'] = (6.4, 4.8)
    
    adata.obs['renamed_samples'] = adata.obs['renamed_samples'].astype('category')
    samples = adata.obs['renamed_samples'].cat.categories
    for sample in samples:
        sample_adata = adata[adata.obs['renamed_samples'] == sample].copy()
        sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].astype('category')
        sample_adata.obs['renamed_samples'] = sample_adata.obs['renamed_samples'].cat.remove_unused_categories()

        # Save into figures/ subfolder
        fig = sc.pl.umap(
            sample_adata,
            color='renamed_samples',
            title=f'Sample: {sample}',
            size=80,
            show=False,
            return_fig=True
        )
        fig.set_size_inches(6.4, 4.8)
        fig.savefig(f'figures/{base}_{sample}_no_leiden.png', bbox_inches='tight', dpi=500)
        plt.close(fig)
    
    # Reset figure size to default (optional)
    plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']

if __name__ == '__main__':
    main()
