#!/usr/bin/env python3
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour- OC Edit
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_mudata', help='Path to MuData object')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')
parser.add_argument('--categorical_covariate_keys', default=None, help='categorical_covariate_keys variables (default: None)')
parser.add_argument('--corrected_counts', default=False, help='Returns corrected counts, aka posterior_predictive_sample() (default: False)')
parser.add_argument('--denoised_data', default=False, help='Returns denoised data, aka get_normalized_expression()  (default: False)')
#parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_mudata = args.path_to_mudata
prefix = args.prefix
batchkey = args.batchkey
#confoundings = args.confoundings
if args.categorical_covariate_keys:
        # Split the string into a list by commas
        categorical_covariate_keys = args.categorical_covariate_keys.split(',')
        #print(categorical_covariate_keys)
else:
    categorical_covariate_keys = None
corrected_counts = args.corrected_counts
denoised_data = args.denoised_data
#latent_key = args.latent_key

# working_dir = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/'
# path_to_mudata = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/export_data/totalvi_igt1_56_20231030_allgenes_mdata.h5mu'
# prefix = 'totalvi_igt1_56_allgenes_20240526_igtsampleregressedout'
# batchkey = 'IGT'
# categorical_covariate_keys = ['sample_id']

# working_dir='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96'
# path_to_mudata='/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/export_data/igt1_96_20241006.h5mu'
# prefix='totalvi_20241006' #prefix = "totalvi_20241008_rmIGTsample"
# batchkey='IGT'
# categorical_covariate_keys='IGTHT'
# corrected_counts=False
# denoised_data=False


print(f"Working Directory: {working_dir}")
print(f"Path to AnnData: {path_to_mudata}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"categorical_covariate_keys: {categorical_covariate_keys}")
print(f"corrected_counts: {corrected_counts}")
print(f"denoised_data: {denoised_data}")
#print(f"Latent Key: {latent_key}")

print("Importing libraries")
import warnings; warnings.simplefilter('ignore')
import scvi
import os
import sys
import scanpy as sc
#import muon as mu
import mudata as mu
import anndata as AnnData
import numpy as np
#import mplscience
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import seaborn as sns
import torch
import fast_matrix_market as fmm
#import pymde #to run MDE

print("Global configurations")
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
scvi.settings.seed = 0  # optional: ensures reproducibility
#pymde.seed(0)
#sns.set_theme()
if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

print("Reading mudata")
os.chdir(working_dir)
#mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
#print(mdata)
#mdata = AnnData.read(working_dir+"/MERGED_RNA_Manual_GSE182509_samples_processed.h5ad") #load this one (results of this block)
mdata = AnnData.read(path_to_mudata)
print(mdata)

#full.mod['RNA'].obs_names = full_RNA.obs_names.astype(str)
#full.mod['RNA'].var_names = genes['x'].astype(str)
#full.mod['ADT'].obs_names = full_ADT.obs_names.astype(str)
#full.mod['ADT'].var_names = proteins.astype(str)

ImmgenT_IGTHT = pd.read_csv(working_dir+"ImmgenT_IGTHT.csv", delimiter = ",")
ImmgenT_IGTHT.index = ImmgenT_IGTHT['Unnamed: 0'].values
ImmgenT_IGTHT = ImmgenT_IGTHT.drop(columns = "Unnamed: 0")


query_IGTHT = pd.read_csv(working_dir+"query_IGTHT.csv", delimiter = ",")
query_IGTHT.index = query_IGTHT['Unnamed: 0'].values
query_IGTHT = query_IGTHT.drop(columns = "Unnamed: 0")
query_IGTHT['Unknown'] = "Unknown"

mdata.obs['IGTHT'] = pd.concat([ImmgenT_IGTHT["IGTHT"],query_IGTHT["IGTHT"]]).astype(str)
mdata.obs['IGT'] = pd.concat([ImmgenT_IGTHT["IGT"],query_IGTHT["IGT"]]).astype(str)
mdata.layers["counts"] = mdata.X.copy()
mdata.obs_names = mdata.obs_names.astype(str)
mdata.var_names = mdata.var_names.astype(str)
mdata.X = mdata.X.copy()
#mdata.var = mdata.var.reset_index(drop=True)
#mdata.mod['RNA'].var = mdata.mod['RNA'].var.reset_index(drop=True)
#mdata.mod['protein'].var = mdata.mod['protein'].var.reset_index(drop=True)
print(mdata)

# print("Filter out genes with no data")
# total_counts_per_genes = np.array(mdata.mod['RNA'].X.sum(axis=0)).flatten()
# low_count_genes = total_counts_per_genes < 1
# print(f"Number of genes with no counts: {low_count_genes.sum()}")
# mdata.mod["RNA"] = mdata.mod["RNA"][:, ~low_count_genes]

print("Creating RNA layer counts")
mdata.layers["counts"] = mdata.X.copy()

#if sparse matrix in protein data: print("Converting  mdata.mod['protein'].X to array")
#mdata.mod['protein'].X = mdata.mod['protein'].X.toarray()

#mdata.update()

print("Batch key list")
print(mdata.obs[batchkey].unique().tolist())
print("Batch key list: any NA?")
print((mdata.obs[batchkey].isna()).any())


print("categorical_covariate_keys")
if categorical_covariate_keys is not None:
    for c in categorical_covariate_keys:
        print(c)
        print(mdata.obs[c].unique().tolist())
        print("any NA?")
        print((mdata.obs[c].isna()).any())
        #print(mdata.mod["RNA"].obs[c].head(10))
        #mdata.mod["RNA"].obs[c] = mdata.mod["RNA"].obs[c].str.replace('.', '', regex=False)
        #print(mdata.mod["RNA"].obs[c].head(10))

print("Setup mu_data")
#scvi.model.TOTALVI.setup_mudata(
#    mdata,
#    rna_layer="counts",
#    categorical_covariate_keys = categorical_covariate_keys, #make sure this is a list!
#    protein_layer=None,
#    batch_key=batchkey,
#    modalities={
#        "rna_layer": "RNA",
#        "protein_layer": "protein",
#        "batch_key": "RNA",
#        "categorical_covariate_keys":"RNA"
#    },
#)

#model = scvi.model.TOTALVI(mdata, n_latent = 30, gene_likelihood = "nb") #extra_decoder_kwargs={'activation_function_bg': 'exp'}


## Use RNA adata instead
#adata = mdata.copy()
#adata
#scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
#model = scvi.model.SCVI(adata, n_latent = 30, gene_likelihood = "nb")
#mdata = adata.copy()

print("Train model")
##uncomment to train model
#model.train()

print("Save model")
#model.save(prefix) #, save_anndata=True)
#model = scvi.model.TOTALVI.load(prefix)
#scvi_model = scvi.model.SCVI.load(prefix, mdata)
# mdata = model.adata
#mdata = model.adata

print("Train new model for SCANVI")
annotation_level2 = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_annotation_table.csv", delimiter = ",", index_col=0)
##annotation_level2.index = annotation_level2['cell_id'].values
##annotation_level2 = annotation_level2.drop(columns = "Unnamed: 0")
annotation_level2 = annotation_level2.loc[annotation_level2.index.intersection(mdata.obs.index)]
##annotation_level2 = annotation_level2[annotation_level2.index.isin([adata.obs.index])]
##mdata.obs["level2"] = pd.concat([annotation_level2["level2"],query_IGTHT["Unknown"]]).astype(str)
mdata.obs['level2'] = mdata.shape[0] * ["Unknown"]
mdata.obs.loc[annotation_level2.index, 'level2'] = annotation_level2['level2'].values

print("Add seed labels")


# Comment if already trained
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
scvi_model = scvi.model.SCVI(mdata, n_latent=30, n_layers=2)
scvi_model.train(50)
scvi_model.save(prefix+"/scvi_model/") #, save_anndata=True)

# Uncomment if model already trained
#scvi_model = scvi.model.SCVI.load(prefix+"/scvi_model/", mdata)
#scvi_model = torch.load(prefix+"/scvi_model/model.pt", weights_only=False)

## Training scanvi model on scvi model
model = scvi.model.SCANVI.from_scvi_model(scvi_model, "Unknown")
model.train()
model.save(prefix+"/scanvi_model/") #, save_anndata=True)

SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "C_scANVI"

mdata.obsm[SCANVI_LATENT_KEY] = model.get_latent_representation(mdata)
mdata.obs[SCANVI_PREDICTIONS_KEY] = model.predict(mdata)
predicted_df = pd.DataFrame(mdata.obs[SCANVI_PREDICTIONS_KEY], index = mdata.obs.index)
predicted_df.to_csv(prefix+"/predicted_celltypes.csv")


##fig, ax = plt.subplots(1, 1)
##model.history["elbo_train"].plot(ax=ax, label="train")
##model.history["elbo_validation"].plot(ax=ax, label="validation")
##ax.set(title="Negative ELBO over training epochs", ylim=(1200, 1400))
##ax.legend()
##fig = plt.gcf() 
##fig.savefig(prefix+"/training_elbo_plot.pdf")

print("Save latent_representation.csv")
latent_representation = scvi_model.get_latent_representation()
TOTALVI_LATENT_KEY = "X_totalVI"
##mdata.mod['RNA'].obsm[TOTALVI_LATENT_KEY] = latent_representation
mdata.obsm[TOTALVI_LATENT_KEY] = latent_representation
##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
latent_df = pd.DataFrame(latent_representation, index = mdata.obs.index)

latent_df.to_csv(prefix+"/latent.csv", index=True)
## latent_df = pd.read_csv(prefix+"/latent.csv", index_col = 0)
## mdata.obsm[TOTALVI_LATENT_KEY] = latent_df
## mdata.obsm[TOTALVI_LATENT_KEY] = mdata.obsm[TOTALVI_LATENT_KEY].values

print("Save umap.csv")
sc.pp.neighbors(mdata, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(mdata, min_dist=0.4)
##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.obs.index)
umap_df.to_csv(prefix+"/umap_python.csv", index=True)

print("Save Anndata")
mdata.write_h5ad(prefix+"/adata_RNA.h5ad")
