import argparse
import scanpy as sc

def main():
    parser = argparse.ArgumentParser(description='Annotate h5ad file with cell type labels')
    parser.add_argument('--input', type=str, required=True, help='Input h5ad file')
    parser.add_argument('--annotations', type=str, required=True, help='Annotation file (csv format: cluster_id,cell_type)')
    parser.add_argument('--output', type=str, required=True, help='Output h5ad file')
    parser.add_argument('--remove_cluster', type=str, required=True, help='Cluster to remove (e.g., "45")')
    
    args = parser.parse_args()
    
    # Read input h5ad
    adata = sc.read_h5ad(args.input)
    
    # Remove specified cluster
    adata = adata[adata.obs["combined_leiden"] != args.remove_cluster].copy()
    
    # Read annotation file
    cluster_to_celltype_dict = {}
    with open(args.annotations, "r") as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                parts = line.split(',')
                if len(parts) == 2:
                    cluster, celltype = parts
                    cluster_to_celltype_dict[cluster] = celltype
    
    cluster_to_celltype_dict = {str(key): value for key, value in cluster_to_celltype_dict.items()}
    
    # Map cluster IDs to cell types
    adata.obs["celltype"] = adata.obs["combined_leiden"].map(cluster_to_celltype_dict)
    
    # Save annotated h5ad
    adata.write_h5ad(args.output)
    
    # UMAP with annotations on top of clusters
    sc.pl.umap(adata, color='celltype', legend_loc="on data", save='_annotationsON.png')
    
    # UMAP with annotations only in legend
    sc.pl.umap(adata, color='celltype', save='_annotations.png')

if __name__ == '__main__':
    main()
