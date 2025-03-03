import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import pandas as pd

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
newObject = "renamed_"+ myObject

combined_adata = sc.read_h5ad(myObject, backed="r")

sample_names = combined_adata.obs['sample']


def rename_samples(row):
    sample_name = row.name 
    sample_type = row["sample"] 

    if sample_type.startswith("TH"):
        print(f"TH Sample: {sample_name}") 
        return "Control"
    
    if sample_type == "Zebra":
        if "-1" in sample_name:
            return "LD"
        elif "-2" in sample_name:
            return "NMDA"

    return sample_type  

combined_adata.obs["renamed_samples"] = combined_adata.obs.apply(rename_samples, axis=1)

print(combined_adata.obs[["sample", "renamed_samples"]].head(10))

combined_adata.write(newObject,compression="gzip")
