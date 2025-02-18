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
        "15": "Proliferating cells", 
        
        "8" : "Rod", 
        "30": "Rod",
        "1" : "Rod", 
        "8" : "Rod",
        "15": "Rod", 
        "17": "Rod",
        "20": "Rod", 
        "36": "Rod", 

        "4" : "MG",
        "25": "MG", 
        "28": "MG", 
        "40": "MG", 
        "48": "MG", 

        "41": "Endothelial",

        "44": "Micorglia", 
        "29": "Micorglia", 
       
        "0" : "Progenitors", 
        "37": "Progenitors", 
        "39": "Progenitors", 
        "45": "Progenitors", 

        "42": "AC", 
        "7" : "AC", 
        "31": "AC", 
        "11": "AC", 
        "10": "AC", 
        "31": "AC", 
        "19": "AC", 
        "24": "AC", 
       
        "23": "BC", 
        "35": "BC",
        "38": "BC",
        "2" : "BC",
        "3" : "BC", 
        "6" : "BC", 
        "9" : "BC", 
        "12": "BC", 
        "13": "BC", 
        "16": "BC", 
        "21": "BC", 
        "22": "BC",
        "26": "BC", 
        "32": "BC",
        "34": "BC", 
        "35": "BC", 


        "5" : "Cone",
        "18": "Cone", 
        "27": "Cone", 
        "43": "Cone",
        "46": "Cone", 

        "14": "RGC", 
        "33": "RGC", 
        "38": "RGC", 


    }
)


combined_adata.obs_names_make_unique()
sc.pl.umap(combined_adata, color='celltype',legend_loc="on data", save="_annotations.png")
combined_adata.write(newObject, compression="gzip")
