import scanpy as sc
import sys
import importlib_metadata
import anndata 

sys.modules['importlib.metadata'] = importlib_metadata


parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")  
newObject = "corrected_" + parts[1]

combined_adata = sc.read(myObject)

batches = combined_adata.obs["sample"].unique()
corrected_batches = []  # List to store corrected batch data

for batch in batches:
    print(f"Processing batch: {batch}")  # Track progress
    adata_batch = combined_adata[combined_adata.obs["sample"] == batch].copy()
    
    # Apply ComBat
    sc.pp.combat(adata_batch, key="sample")
    
    # Append corrected batch to list
    corrected_batches.append(adata_batch)

corrected_adata = anndata.concat(corrected_batches, join="outer", merge="same")

corrected_adata.write(newObject)


