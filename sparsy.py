import scanpy as sc
import argparse
import scipy.sparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Downsample and filter an h5ad file.")
parser.add_argument('--input', '-i', required=True, help="Path to input .h5ad file")
parser.add_argument('--output', '-o', required=True, help="Path to output .h5ad file")
args = parser.parse_args()

# Read input file
adata = sc.read_h5ad(args.input)

# Keep only highly variable genes and downsample cells
adata = adata[:, adata.var['highly_variable']]
adata = adata[::5, :]  # downsample every 5th cell

# Ensure sparse format
adata.X = scipy.sparse.csr_matrix(adata.X)

# Write output file
sc.write(args.output, adata)

