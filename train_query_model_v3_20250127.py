#!/usr/bin/env python3

# python3 $1/train_query_model_v3.py [query_data_path] [protein_path] [integrated_data_path] [subgroup]

import pandas as pd
import mudata as mu
import anndata as AnnData
import scanpy as sc
from scipy.sparse import csc_matrix
import scvi
import numpy as np
import scanpy
import matplotlib.pyplot as plt
import mudata
from scipy.io import mmread
import fast_matrix_market as fmm
from itertools import chain
import torch
import sys
import time

prefix = "Trial_TOTALVI_PostQC"
categorical_covariate_keys = ['IGTHT']
batchkey = 'IGT'


## read in data
#query = mudata.read("IGT97_104_mudata.h5mu")
#ref = mudata.read("igt1_96_20250109.h5mu")
#query_metadata = query.obs
#ref_metadata = ref.obs
#all_rna = AnnData.concat([ref.mod['RNA'], query.mod['RNA']])

## Protein Data
#protein_data = np.zeros((len(query.obs_names), 180))
#protein_data = protein_data.astype(int)
#protein_df = pd.DataFrame(protein_data)
#protein_df.index = query.obs_names
#protein_df.columns = ref.mod['protein'].var_names
#protein_query_adata = AnnData.AnnData(protein_df)
#all_protein = AnnData.concat([ref.mod['protein'], protein_query_adata])
#all_mudata = mu.MuData({"RNA": all_rna, "protein": all_protein})


#all_metadata = pd.concat([ref_metadata, query_metadata])
#all_mudata.obs = all_metadata
#all_mudata.mod['RNA'].obs = all_mudata.obs

#all_metadata.to_csv("all_metadata.csv", index = 1)
#all_mudata.obs.drop(columns=["is_spleen_standard"], inplace=True)
#all_mudata.write_h5mu("all_mudata.h5mu")

#all_mudata.mod['RNA'].layers["counts"] = all_mudata.mod['RNA'].X.copy()
#all_mudata.update()

## Read in already created mudata object
all_mudata = mu.read("all_data_PostQC.h5mu")


## setup mudata
print("Setup mu_data")
scvi.model.TOTALVI.setup_mudata(
    all_mudata,
    rna_layer="counts",
    categorical_covariate_keys = categorical_covariate_keys, #make sure this is a list!
    protein_layer=None,
    batch_key=batchkey,
    modalities={
        "rna_layer": "RNA",
        "protein_layer": "protein",
        "batch_key": "RNA",
        "categorical_covariate_keys":"RNA"
    },
)
## Initialize model
model = scvi.model.TOTALVI(all_mudata, n_latent = 30, gene_likelihood = "nb")

## Train model
print("Train model")
start_time = time.time()
print(f"Start time: {start_time:.2f} seconds")
model.train()

end_time = time.time()
print(f"End time: {end_time:.2f} seconds")

elapsed_time = end_time - start_time
print(f"Total elapsed time: {elapsed_time:.2f} seconds")


## Save model
print("Save model")
model.save(prefix+"/model/", save_anndata=True)

## get latent representation
print("Save latent_representation.csv")
latent_representation = model.get_latent_representation()
TOTALVI_LATENT_KEY = "X_totalVI"
all_mudata.obsm[TOTALVI_LATENT_KEY] = latent_representation
latent_df = pd.DataFrame(latent_representation, index = all_mudata.mod['RNA'].obs.index)
latent_df.to_csv(prefix+"/latent.csv", index=True)

print("Save umap.csv")
sc.pp.neighbors(all_mudata, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(all_mudata, min_dist=0.4)
umap_df = pd.DataFrame(all_mudata.obsm['X_umap'], index = all_mudata.mod['RNA'].obs.index)
umap_df.to_csv(prefix+"/umap_python.csv", index=True)

sc.pl.umap(
    all_mudata,
    color=["IGT"],
    frameon=False,
    ncols=1,
    title="Query & reference",
    save = prefix+"/all_mudata_umap_IGT.png"
)

sc.pl.umap(
    all_mudata,
    color=["query"],
    frameon=False,
    ncols=1,
    title="Query & reference",
    save = prefix+"/all_mudata_umap_query.png"
)

