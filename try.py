import scanpy as sc
import sys 
import importlib_metadata

sys.modules['importlib.metadata'] = importlib_metadata

samples = {"Zebra":"Zebra_filtered_feature_bc_matrix.h5",
        "TH115": "TH115_filtered_feature_bc_matrix.h5",
        "TH44":"TH44_filtered_feature_bc_matrix.h5",
        "TH54":"TH54_filtered_feature_bc_matrix.h5",
        "TH55":"TH55_filtered_feature_bc_matrix.h5",
        "TH56":"TH56_filtered_feature_bc_matrix.h5",
        "TH57":"TH57_filtered_feature_bc_matrix.h5",
        "TH71":"TH71_filtered_feature_bc_matrix.h5" }


adatas = {}

for sample_id, filename in samples.items():
    print(f"Reading {filename}...")  # Optional: Print status
    adata = sc.read_10x_h5(filename)  # Read the file
    adata.var_names_make_unique() 
    adata.obs['sample'] = sample_id 
    adatas[sample_id] = adata 

combined_adata = sc.concat(adatas.values(), label='sample', keys=adatas.keys())

combined_adata.var["mt"] = adata.var_names.str.startswith("mt-")

sc.pl.violin(
    combined_adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)


sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)

sc.pp.highly_variable_genes(combined_adata, flavor='seurat', n_top_genes=2000)

sc.pp.scale(combined_adata, max_value=10)

sc.tl.pca(combined_adata, svd_solver='arpack')

sc.pp.neighbors(combined_adata)

sc.tl.umap(combined_adata)

sc.pl.umap(combined_adata, color='sample', size=2, save='umap_combined.png')





