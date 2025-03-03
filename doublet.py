import scrublet as scr
import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import scipy.sparse as sp


sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
args = parser.parse_args()

myObject =  args.myObject
parts = myObject.split("_")

newObject = "doubletRemoved_" + parts[1]
parts = myObject.split("_")

fname = parts[1].split(".") 

fname = fname[0] 

combined_adata = sc.read(myObject, backed="r")

counts_matrix = combined_adata.X[:] 

scrub = scr.Scrublet(counts_matrix)

doublet_scores, predicted_doublets = scrub.scrub_doublets()

combined_adata.obs['predicted_doublets'] = predicted_doublets

figure_name = fname + "predictedDoublets.png"
sc.pl.violin(combined_adata, 'predicted_doublets', jitter=True, boxplot=True, save = figure_name)

combined_adata = combined_adata[~combined_adata.obs['predicted_doublets'], :]

combined_adata.write(newObject, compression="gzip")

