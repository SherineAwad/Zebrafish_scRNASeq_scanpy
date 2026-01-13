#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--annotations', required=True)
    args = parser.parse_args()

    # Read
    adata = sc.read_h5ad(args.input)

    # Read annotations file
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

    # Create new column
    adata.obs['combined_leiden'] = leiden_str.map(lambda x: mapping.get(x, x))

    # Save
    adata.write_h5ad(args.output)

    # Plots
    base = args.output.replace('.h5ad', '')

    # 1. UMAP colored by sample
    sc.pl.umap(adata, color='renamed_samples', save=f'{base}_colored_by_sample.png')

    # 2. UMAP colored by combined_leiden (all samples together)
    sc.pl.umap(adata, color='combined_leiden', save=f'{base}_combined_leiden.png')

    # 3. Separate UMAP for each sample with combined_leiden
    for sample in adata.obs['renamed_samples'].unique():
        mask = adata.obs['renamed_samples'] == sample
        adata_sample = adata[mask].copy()
        sc.pl.umap(adata_sample, color='combined_leiden', title=f'Sample: {sample}', save=f'{base}_{sample}.png')
    
    # 4. Separate UMAP for each sample without leiden (just points)
    for sample in adata.obs['renamed_samples'].unique():
        mask = adata.obs['renamed_samples'] == sample
        adata_sample = adata[mask].copy()
        sc.pl.umap(adata_sample, color=None, title=f'Sample: {sample}', save=f'{base}_{sample}_no_leiden.png')

if __name__ == '__main__':
    main()
