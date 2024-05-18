#!/project/zemmour/david/envs/scvi_20240315/bin/python
"""This script run SCVI from a Anndata object specifying batch_key and categorical_covariate_keys if needed"""
#author: David Zemmour
#date: 05/17/2024
#run_scvi.py [cwd] [path to h5ad] [prefix] [batchkey] [confounding] [name of latent space]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run SCVI from a Anndata object specifying batch_key and categorical_covariate_keys if needed")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_adata', help='Path to AnnData object')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', help='Batch key for analysis')
parser.add_argument('--confoundings', default=None, help='Confounding variables (default: None)')
parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_adata = args.path_to_adata
prefix = args.prefix
batchkey = args.batchkey
#confoundings = args.confoundings
if args.confoundings:
        # Split the string into a list by commas
        confoundings = args.confoundings.split(',')
        print(confoundings)
else:
    confoundings = None
latent_key = args.latent_key

print(f"Working Directory: {working_dir}")
print(f"Path to AnnData: {path_to_adata}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"Confounding Variable: {confoundings}")
print(f"Latent Key: {latent_key}")

print("Importing libraries")
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

print("Read adata")
os.chdir(working_dir)
adata = sc.read_h5ad(path_to_adata)

#Setup adata
print("Setup anndata")

scvi.model.SCVI.setup_anndata(
    adata,
    categorical_covariate_keys = confoundings, #make sure this is a list!
    batch_key=batchkey
)

model = scvi.model.SCVI(adata, n_latent = 30, n_hidden =128, n_layers = 2, gene_likelihood = "nb")

print("Train model")
model.train()

print("Save model")
model.save(prefix, save_anndata=True)

# fig, ax = plt.subplots(1, 1)
# model.history["elbo_train"].plot(ax=ax, label="train")
# model.history["elbo_validation"].plot(ax=ax, label="validation")
# ax.set(title="Negative ELBO over training epochs", ylim=(1200, 1400))
# ax.legend()
# fig.savefig(prefix+"/training_elbo_plot.pdf")

latent_representation = model.get_latent_representation()
adata.obsm[latent_key] = latent_representation
latent_df = pd.DataFrame(latent_representation, index = adata.obs.index)

latent_df.to_csv(prefix+"/latent_space.csv", index=True)
adata.write(prefix+"/adata.h5ad")
