import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import anndata
import scipy.sparse as sp
import numpy as np



sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")  
newObject = "annotated_" + parts[1]

combined_adata = sc.read_h5ad(myObject, backed="r")

combined_adata.obs["celltype"] = combined_adata.obs["leiden"].map(
    {

        "0" : "MG", 
        "4" : "MG",
        "25": "MG",
    
        "15": "Progenitors", 

        "8": "Post Mitotic Rod Precursors",

        "17": "Rod", 
        "20": "Rod", 
        
        "1": "Cones", 
        "5": "Cones", 
        "18": "Cones", 

        "2" : "BC", 
        "3" : "BC",
        "6": "BC",
        "9": "BC",
        "12": "BC",
        "13": "BC",
        "16": "BC",
        "21": "BC",
        "22": "BC",
        "23": "BC",
        "26": "BC",
        "30": "BC",
        "32": "BC",
        "34": "BC",
        "35": "BC",
        "38": "BC",
        
        "7": "AC", 
        "24":"AC",
        "27":"AC",
        "33":"AC",
        "37":"AC",
        "39":"AC",
        "42":"AC",

        "14": "RGC", 
        "29": "Microglia", 

        "41": "Endothelial Cells", 

        "28": "RPE", 
        "40": "RPE", 
     
        "11": "HC", 
        "19": "HC", 
        "10": "HC", 
        "31": "HC" 

    }
)


combined_adata.obs_names_make_unique()
sc.pl.umap(combined_adata, color='celltype',legend_loc="on data", save="_annotations.png")
combined_adata.write(newObject, compression="gzip")
