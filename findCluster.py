import scanpy as sc
import copy 
import pandas as pd 

combined_adata = sc.read_h5ad("annotated_Zebrafishes.h5ad")

combined_adata.uns['log1p'] = {'base': 1.0}  # Assuming log1p with base 1 (natural log)

sc.tl.rank_genes_groups(combined_adata, groupby='leiden', groups=['53'], reference='rest', method='wilcoxon')
cluster_53_results = copy.deepcopy(combined_adata.uns['rank_genes_groups'])

top_genes_53 = pd.DataFrame({
    'gene': combined_adata.uns['rank_genes_groups']['names']['53'],
    'score': combined_adata.uns['rank_genes_groups']['scores']['53'],
    'pval': combined_adata.uns['rank_genes_groups']['pvals']['53']
})

print(top_genes_53.head(10))

# Get the top 10 genes
top_genes_list = top_genes_53['gene'].head(10).tolist()

# Plot UMAP for top genes and save the plot
sc.pl.umap(combined_adata, color=top_genes_list, size=50, cmap='viridis', save='_top_genes_53_umap.png')



sc.tl.rank_genes_groups(combined_adata, groupby='leiden', groups=['55'], reference='rest', method='wilcoxon')
cluster_55_results = copy.deepcopy(combined_adata.uns['rank_genes_groups'])

top_genes_55 = pd.DataFrame({
    'gene': combined_adata.uns['rank_genes_groups']['names']['55'],
    'score': combined_adata.uns['rank_genes_groups']['scores']['55'],
    'pval': combined_adata.uns['rank_genes_groups']['pvals']['55']
})

print(top_genes_55.head(10))

# Get the top 10 genes
top_genes_list = top_genes_55['gene'].head(10).tolist()

# Plot UMAP for top genes and save the plot
sc.pl.umap(combined_adata, color=top_genes_list, size=50, cmap='viridis',  save='_top_genes_55_umap.png')


