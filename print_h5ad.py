#!/usr/bin/env python3

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse

def summarize_matrix(matrix, name, gene_names=None, cell_names=None):
    """Helper function to summarize a matrix with statistics"""
    print(f"\n=== {name} ===")
    
    if matrix is None:
        print("None")
        return
    
    # Convert to dense if sparse for certain operations
    if issparse(matrix):
        dense_matrix = matrix.toarray()
        matrix_type = "sparse"
    else:
        dense_matrix = matrix
        matrix_type = "dense"
    
    print(f"Shape: {matrix.shape} (cells: {matrix.shape[0]}, genes: {matrix.shape[1]})")
    print(f"Type: {matrix_type}")
    print(f"Data type: {matrix.dtype}")
    
    # Basic statistics
    print(f"Min: {dense_matrix.min():.6f}")
    print(f"Max: {dense_matrix.max():.6f}")
    print(f"Mean: {dense_matrix.mean():.6f}")
    print(f"Std: {dense_matrix.std():.6f}")
    
    # Count zeros and non-zeros
    if issparse(matrix):
        non_zero_count = matrix.nnz
        zero_count = matrix.shape[0] * matrix.shape[1] - non_zero_count
        sparsity = zero_count / (matrix.shape[0] * matrix.shape[1])
    else:
        non_zero_count = np.count_nonzero(dense_matrix)
        zero_count = dense_matrix.size - non_zero_count
        sparsity = zero_count / dense_matrix.size
    
    print(f"Non-zero values: {non_zero_count}")
    print(f"Zero values: {zero_count}")
    print(f"Sparsity: {sparsity:.4f} ({sparsity*100:.2f}%)")
    
    # Show top genes and cells
    print(f"\nTop 5 genes (first 5 rows):")
    if gene_names is not None and len(gene_names) > 0:
        for i in range(min(5, len(gene_names))):
            gene_data = dense_matrix[:, i] if dense_matrix.shape[1] > i else []
            print(f"  {gene_names[i]}: min={gene_data.min():.4f}, max={gene_data.max():.4f}, mean={gene_data.mean():.4f}")
    
    print(f"\nTop 5 cells (first 5 columns):")
    if cell_names is not None and len(cell_names) > 0:
        for i in range(min(5, len(cell_names))):
            cell_data = dense_matrix[i, :] if dense_matrix.shape[0] > i else []
            print(f"  {cell_names[i]}: min={cell_data.min():.4f}, max={cell_data.max():.4f}, mean={cell_data.mean():.4f}")

def main():
    # ------------------------
    # Parse command-line args
    # ------------------------
    parser = argparse.ArgumentParser(
        description="Print detailed information about .h5ad file including categorical observations and expression matrices"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the annotated .h5ad file"
    )
    args = parser.parse_args()

    # ------------------------
    # Load the AnnData object
    # ------------------------
    print(f"Loading AnnData object from: {args.input_file}")
    adata = sc.read_h5ad(args.input_file)
    
    # Basic AnnData information
    print(f"\n{'='*50}")
    print("BASIC AnnData INFORMATION")
    print(f"{'='*50}")
    print(f"Overall shape: {adata.shape} (cells: {adata.n_obs}, genes: {adata.n_vars})")
    print(f"Observation names (cells): {list(adata.obs_names[:5])}..." if len(adata.obs_names) > 5 else f"Observation names: {list(adata.obs_names)}")
    print(f"Variable names (genes): {list(adata.var_names[:5])}..." if len(adata.var_names) > 5 else f"Variable names: {list(adata.var_names)}")

    # ------------------------
    # Summarize adata.X
    # ------------------------
    summarize_matrix(
        adata.X, 
        "adata.X (Primary Expression Matrix)",
        gene_names=list(adata.var_names),
        cell_names=list(adata.obs_names)
    )

    # ------------------------
    # Summarize adata.raw
    # ------------------------
    if adata.raw is not None:
        raw_adata = adata.raw
        print(f"\nadata.raw type: {type(raw_adata)}")
        
        if hasattr(raw_adata, 'X') and raw_adata.X is not None:
            summarize_matrix(
                raw_adata.X,
                "adata.raw.X (Raw Expression Matrix)",
                gene_names=list(raw_adata.var_names) if hasattr(raw_adata, 'var_names') else None,
                cell_names=list(raw_adata.obs_names) if hasattr(raw_adata, 'obs_names') else None
            )
            
            # Compare adata.X and adata.raw.X if both exist
            if adata.X is not None and raw_adata.X is not None:
                print(f"\n{'='*50}")
                print("COMPARISON: adata.X vs adata.raw.X")
                print(f"{'='*50}")
                print(f"Same shape: {adata.X.shape == raw_adata.X.shape}")
                if adata.X.shape == raw_adata.X.shape:
                    if issparse(adata.X) and issparse(raw_adata.X):
                        diff = (adata.X != raw_adata.X).sum()
                    else:
                        diff = np.sum(adata.X != raw_adata.X)
                    print(f"Number of different elements: {diff}")
                    print(f"Identical matrices: {diff == 0}")
        else:
            print("\nadata.raw exists but adata.raw.X is None")
    else:
        print("\nNo adata.raw found")

    # ------------------------
    # Identify categorical columns
    # ------------------------
    print(f"\n{'='*50}")
    print("CATEGORICAL OBSERVATIONS")
    print(f"{'='*50}")
    
    categorical_cols = [
        col for col in adata.obs.columns
        if adata.obs[col].dtype.name == 'category' or adata.obs[col].dtype == object
    ]

    if not categorical_cols:
        print("No categorical columns found in this dataset.")
    else:
        print(f"Found {len(categorical_cols)} categorical columns:\n")
        for col in categorical_cols:
            unique_values = adata.obs[col].unique()
            print(f"Column: {col} ({len(unique_values)} unique values)")
            print(f"Values: {list(unique_values)}")
            print("---")

    # ------------------------
    # Additional metadata information
    # ------------------------
    print(f"\n{'='*50}")
    print("ADDITIONAL METADATA")
    print(f"{'='*50}")
    print(f"adata.obs columns: {list(adata.obs.columns)}")
    print(f"adata.var columns: {list(adata.var.columns)}")
    
    # Check for layers
    if adata.layers:
        print(f"\nLayers available: {list(adata.layers.keys())}")
        for layer_name, layer_matrix in adata.layers.items():
            summarize_matrix(
                layer_matrix,
                f"Layer: {layer_name}",
                gene_names=list(adata.var_names),
                cell_names=list(adata.obs_names)
            )
    else:
        print("\nNo additional layers found.")

    print(f"\n{'='*50}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
