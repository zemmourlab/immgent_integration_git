#!/project/zemmour/david/envs/scvi_20240315/bin/python
"""This script run TOTALVI denoised data"""
#author: David Zemmour
#date: 04/12/2024
#run_totalvi_denoising_v2.py [cwd] [prefix]

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
prefix = sys.argv[2] #"totalvi_igt1_56_allgenes_Treg_20240327_organregressedout"
# batchkey = sys.argv[4]
# confounding1 = sys.argv[5]
# totalvi_latent_key = sys.argv[6]

# working_dir = "/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56"
# prefix = "totalvi_igt1_56_20231030_allgenes_organregressed"

print("working_dir:"+working_dir)
print("prefix:"+prefix)

model = scvi.model.TOTALVI.load(prefix)
mdata = model.adata

print("Calculate denoised protein data")
denoised = model.get_normalized_expression(gene_list = []) 

mdata.mod['protein'].layers["counts"] = mdata.mod["protein"].X.copy()
mdata.mod['protein'].layers["protein_denoised"] = denoised[1]
#mdata.mod['RNA'].layers["rna_denoised"] = denoised[0]

print("Save denoised data")
mdata.write(prefix+"/adata_withdenoised.h5mu")

#export mu data in separate AnnData for R
print("Export data rna_data.h5ad and protein_data.h5ad in AnnData for R")
mdata.mod['RNA'].write_h5ad(prefix+"/rna_data.h5ad")
mdata.mod['protein'].write_h5ad(prefix+"/protein_data.h5ad")

print("Done")
