import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")  
newObject = "annotated_" + parts[1]

combined_adata = sc.read(myObject)

combined_adata.obs["celltype"] = combined_adata.obs["leiden"].map(
    {
        "13": "HC",
        "26": "HC", 
        "29": "Microglia",
        "4": "RGC", 
        "0" : "MG", 
        "1" : "Cones", 
        "18": "Cones", 
        "20": "Cones",  
        "17": "Rod", 
        "12": "BC", 
        "22": "BC", 
        "6" : "BC", 
        "3" : "BC",
        "2" : "BC", 
        "16": "BC",
        "9" : "BC", 
        "21": "BC",
        "36": "BC",  
        "34": "BC",
        "32": "BC",
        "24": "AC",
        "7" : "AC",
        "14": "AC", 
        "33": "AC",  
        "42": "AC",
        "41": "Endothelial", 
        "37": "Progenitors", 
        "39": "Progenitors", 
        "45": "Progenitors", 
        "0": "Progenitors", 
        "8" :"Progenitors", 
        "15" : "Progenitors"

    }
)

combined_adata.obs_names_make_unique()
sc.pl.umap(combined_adata, color='celltype',legend_loc="on data", save="_annotations.png")
combined_adata.write(newObject)


