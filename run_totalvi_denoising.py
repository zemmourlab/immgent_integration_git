# David Zemmour
# python 3.9
# usage: python run_totalvi_denoising.py [working directory] [prefix, also serves as path to the model]
#python run_totalvi_denoising.py /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration totalvi_20230508

print("Importing libraries")
import os
import sys
import tempfile

try:
    import scanpy as sc
except Exception as e:
    print("scanpy error usually not a problem undefined symbol: cblas_cdotc_sub")

import scvi
import torch
import scanpy as sc
import anndata as ad

os.chdir(sys.argv[1])

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

#sc.set_figure_params(figsize=(4, 4))
torch.set_float32_matmul_precision("high")

prefix = sys.argv[2]
#prefix = "totalvi_igt1_48_20230706_allgenes"
#prefix = "totalvi_20230508" #spleen standard only, 2000 genes as a test

print("Loading the model")
model = scvi.model.TOTALVI.load(prefix+"/")

print("Denoising")
denoised = model.get_normalized_expression() 
#denoised = model.get_normalized_expression(gene_list = ["Cd8a", "Cd4"], protein_list = ['CD3', 'CD4']) #didn't work: TypeError: tuple indices must be integers or slices, not list

print("Adding denoised data to anndata object")
#add denoised RNA to layer
model.adata.layers["X.denoised"] = denoised[0]
#add denoised protein in obsm["protein_denoised"]
model.adata.obsm["protein_denoised"] = denoised[1]

print("AWriting new anndata object prefix_denoised.h5ad")
model.adata.write(prefix+"denoised.h5ad", compression="gzip")

