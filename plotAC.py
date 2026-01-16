#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

def get_filename_prefix(filepath):
    """
    Extract the part before the first underscore from filename.
    Example: 
      '/path/to/AC_sample_data.h5ad' -> 'AC'
      'AC_experiment1.h5ad' -> 'AC'
      'AC_only.h5ad' -> 'AC'
    """
    filename = os.path.basename(filepath)
    # Remove extension
    name_without_ext = filename.rsplit('.', 1)[0]
    # Get part before first underscore
    if '_' in name_without_ext:
        return name_without_ext.split('_')[0]
    else:
        return name_without_ext  # Return whole name if no underscore

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    args = parser.parse_args()

    # Get AC prefix from input filename (part before first underscore)
    base = get_filename_prefix(args.input)
    print(f"Using prefix: {base}")
    
    # Read input AnnData
    adata = sc.read_h5ad(args.input)

    os.makedirs('figures', exist_ok=True)

    # 1. UMAP colored by sample
    sc.pl.umap(adata, color='renamed_samples', size=80, save=f'{base}_colored_by_sample.png', show=False)

    # 2. UMAP colored by combined_leiden
    sc.pl.umap(adata, color='combined_leiden', size=80, legend_loc='on data', save=f'{base}_combined_leiden.png', show=False)

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
            size=80,legend_loc='on data',
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
