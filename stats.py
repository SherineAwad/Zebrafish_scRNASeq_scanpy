#!/usr/bin/env python3
"""
Calculate cell ratios from h5ad file and save stacked bar plot.
"""

import argparse
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

def main():
    parser = argparse.ArgumentParser(description='Calculate cell ratios from h5ad file')
    parser.add_argument('--input', '-i', required=True, help='Input h5ad file')
    parser.add_argument('--column', '-c', required=True, help='Cell type column name in adata.obs')
    args = parser.parse_args()
    
    # Load data
    adata = sc.read_h5ad(args.input)
    
    # Get prefix from input filename
    prefix = os.path.splitext(os.path.basename(args.input))[0]
    
    # Calculate cell counts per sample and cluster
    ct_counts = pd.crosstab(adata.obs['renamed_samples'], adata.obs[args.column])
    
    # Calculate ratios
    ratios_df = ct_counts.div(ct_counts.sum(axis=1), axis=0)
    
    # Create a combined DataFrame for CSV output
    with open(f"{prefix}_cell_ratios.csv", 'w') as f:
        # Write ratios
        f.write("Cell Ratios (Fractions)\n")
        ratios_df.to_csv(f)
        f.write("\n\n")
        
        # Write counts
        f.write("Cell Counts\n")
        ct_counts.to_csv(f)
        f.write("\n\n")
        
        # Write total cells per sample
        f.write("Total Cells per Sample\n")
        total_per_sample = ct_counts.sum(axis=1)
        total_per_sample.to_csv(f, header=False)
        f.write("\n\n")
        
        # Write total cells per cluster
        f.write("Total Cells per Cluster\n")
        total_per_cluster = ct_counts.sum(axis=0)
        total_per_cluster.to_csv(f, header=False)
    
    print(f"Saved full data to: {prefix}_cell_ratios.csv")
    
    # Print summary to console
    print("\n=== SUMMARY ===")
    print(f"\nTotal cells per sample:")
    total_per_sample = ct_counts.sum(axis=1)
    for sample, count in total_per_sample.items():
        print(f"  {sample}: {count} cells")
    
    print(f"\nTotal cells per cluster:")
    total_per_cluster = ct_counts.sum(axis=0)
    for cluster, count in total_per_cluster.items():
        print(f"  {cluster}: {count} cells")
    
    print(f"\nCell Ratios (first 10 clusters):")
    print(ratios_df.iloc[:, :10])
    print(f"\nFull data saved to {prefix}_cell_ratios.csv")
    
    # Create colors for plot
    colors = plt.cm.tab20.colors
    
    # Plot
    ax = ratios_df.plot(kind='bar', stacked=True, figsize=(12, 6), color=colors)
    plt.ylabel("Fraction of cells")
    plt.xlabel("Sample")
    plt.title("Cell Type Contribution per Sample")
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(f"{prefix}_cell_ratios.png", dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"\nSaved plot to: {prefix}_cell_ratios.png")

if __name__ == "__main__":
    main()
