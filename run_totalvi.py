#!/project/zemmour/david/envs/scvi_20241008/bin/python
"""This script run TOTALVI from a muon object specifying batch_key and categorical_covariate_keys if needed"""
#author: David Zemmour
#date: 04/10/2024
#run_totalvi.py [cwd] [path to mdata] [prefix] [batchkey] [confounding] [name of latent space]

import warnings; warnings.simplefilter('ignore')

import scvi
import os
import sys
import scanpy as sc
import muon as mu

import numpy as np
import mplscience
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch

# Global configurations
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
sc.set_figure_params(figsize=(6, 6), frameon=False)
#sns.set_theme()

if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

#%config InlineBackend.print_figure_kwargs={"facecolor": "w"}
#%config InlineBackend.figure_format="retina"

print("Arguments")
working_dir = sys.argv[1]
path_to_mdata = sys.argv[2]
prefix = sys.argv[3] #"totalvi_igt1_56_allgenes_Treg_20240327_organregressedout"
batchkey = sys.argv[4]
confounding1 = sys.argv[5]
totalvi_latent_key = sys.argv[6]

# working_dir = "/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56"
# path_to_mdata = "/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/export_data/totalvi_igt1_56_20231030_allgenes_mdata.h5mu"
# prefix = "totalvi_igt1_56_20231030_allgenes_organregressed"
# batchkey = "IGT"
# confounding1 = "organ_simplified"
# totalvi_latent_key = "X_totalvi_rmorgan"

print("working_dir:"+working_dir)
print("path_to_mdata:"+path_to_mdata)
print("prefix:"+prefix)
print("batchkey:"+batchkey)
print("confounding1:"+confounding1)
print("totalvi_latent_key:"+totalvi_latent_key)

#Read mdata
print("Read mdata")
os.chdir(working_dir)
mdata = mu.read(path_to_mdata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
print(mdata)

print(mdata.mod["RNA"].obs[batchkey].unique().tolist())
print(mdata.mod["RNA"].obs[confounding1].unique().tolist())

#Setup mu_data
print("Setup mu_data")
scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    categorical_covariate_keys = [confounding1], #make sure this is a list!
    protein_layer=None,
    batch_key=batchkey,
    modalities={
        "rna_layer": "RNA",
        "protein_layer": "protein",
        "batch_key": "RNA",
        "categorical_covariate_keys":"RNA"
    },
)

model = scvi.model.TOTALVI(mdata, n_latent = 30, gene_likelihood = "nb")

print("Train model")
model.train()

print("Save model")
model.save(prefix, save_anndata=True)

fig, ax = plt.subplots(1, 1)
model.history["elbo_train"].plot(ax=ax, label="train")
model.history["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs", ylim=(1200, 1400))
ax.legend()
fig.savefig(prefix+"/training_elbo_plot.pdf")

latent_representation = model.get_latent_representation()
mdata.mod['RNA'].obsm[totalvi_latent_key] = latent_representation
latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)

print("Save latent data in csv file")
latent_df.to_csv(prefix +'/'+  prefix + '_' + totalvi_latent_key + ".csv", index=True)

print("Save adata.h5mu")
mdata.write(prefix+"/adata.h5mu")

# print("Calculate denoised data")
# denoised = model.get_normalized_expression() 
# 
# mdata.mod['protein'].layers["counts"] = mdata.mod["protein"].X.copy()
# mdata.mod['protein'].layers["protein_denoised"] = denoised[1]
# mdata.mod['RNA'].layers["rna_denoised"] = denoised[0]
# 
# print("Save denoised data")
# mdata.write(prefix+"/adata_withdenoised.h5mu")
# 
# #export mu data in separate AnnData for R
# print("Export data rna_data.h5ad and protein_data.h5ad in AnnData for R")
# mdata.mod['RNA'].write_h5ad(prefix+"/rna_data.h5ad")
# mdata.mod['protein'].write_h5ad(prefix+"/protein_data.h5ad")
# 
# print("Done")
