#!/usr/bin/env python
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: Odhran Casey
#date: 2/28/2025
#run_SCANVI_trial.py [cwd] [path to adata .h5ad] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

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
import os
import sys
import scarches as sca
import scanpy as sc
#import mudata as mu
import anndata as AnnData
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch

sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

print("Global configurations")
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
#sca.settings.seed = 0  # optional: ensures reproducibility
torch.manual_seed(0)
#pymde.seed(0)
#sns.set_theme()
if torch.cuda.is_available():
    print("CUDA is available")
    print("Using device:", torch.cuda.get_device_name())
    torch.set_float32_matmul_precision("high")

print("Reading mudata")
os.chdir(working_dir)
adata = AnnData.read(working_dir+"MERGED_RNA_gdT_query_ImmgenT.h5ad") #load this one (results of this block)
cells = adata.obs_names.astype(str)

ImmgenT_IGTHT = pd.read_csv(working_dir+"ImmgenT_IGTHT.csv", delimiter = ",")
ImmgenT_IGTHT.index = ImmgenT_IGTHT['Unnamed: 0'].values
ImmgenT_IGTHT = ImmgenT_IGTHT.drop(columns = "Unnamed: 0")
ImmgenT_annotation =  pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/annotation_table_20250109.csv", delimiter = ",")
ImmgenT_annotation.index = ImmgenT_annotation['cell_id'].values
ImmgenT_annotation = ImmgenT_annotation.drop(columns = "cell_id")
#ImmgenT_annotation = ImmgenT_annotation[ImmgenT_annotation["level1"] == "gdT"]
ImmgenT_level2 = ImmgenT_annotation["level2"]
ImmgenT_annotation

query_IGTHT = pd.read_csv(working_dir+"query_IGTHT.csv", delimiter = ",", index_col=0)
#query_IGTHT.index = query_IGTHT['Unnamed: 0'].values
#query_IGTHT = query_IGTHT.drop(columns = "Unnamed: 0")
query_metadata =  pd.read_csv(working_dir+"/cell_metadata.csv", delimiter = ",", index_col = 0)
#query_metadata.index = query_metadata['Unnamed: 0'].values
#query_metadata = query_metadata.drop(columns = "Unnamed: 0")
query_annotation = query_metadata["celltypes"]
query_annotation

adata.obs['IGTHT'] = pd.concat([ImmgenT_IGTHT["IGTHT"],query_IGTHT["IGTHT"]]).astype(str)
adata.obs['IGT'] = pd.concat([ImmgenT_IGTHT["IGT"],query_IGTHT["IGT"]]).astype(str)
adata.layers["counts"] = adata.X.copy()
adata.obs_names = adata.obs_names.astype(str)
adata.var_names = adata.var_names.astype(str)
adata.X = adata.X.copy()
print(adata)

print("Creating layer counts")
adata.layers["counts"] = adata.X.copy()
#adata.update()

print("Batch key list")
print(adata.obs[batchkey].unique().tolist())
print("Batch key list: any NA?")
print((adata.obs[batchkey].isna()).any())


print("categorical_covariate_keys")
if categorical_covariate_keys is not None:
    for c in categorical_covariate_keys:
        print(c)
        print(adata.obs[c].unique().tolist())
        print("any NA?")
        print((adata.obs[c].isna()).any())
        #print(adata.mod["RNA"].obs[c].head(10))
        #adata.obs[c] = adata.obs[c].str.replace('.', '', regex=False)
        #print(adata.obs[c].head(10))


# SCANVI
# Separate reference and query adata
source_adata =  adata[ImmgenT_annotation.index]
print(source_adata)
target_adata = adata[query_metadata.index]
print(target_adata)
source_adata.obs["level2"] = ImmgenT_annotation["level2"]
source_adata.obs["reference/query"] = "ImmgenT_gdT"
target_adata.obs["celltypes"] = query_metadata["celltypes"]
target_adata.obs["reference/query"] = "query"

# setup and train reference adata
cell_type_key = "level2"
sca.models.SCVI.setup_anndata(source_adata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys, labels_key=cell_type_key)
vae = sca.models.SCVI(source_adata, n_latent = 30, gene_likelihood = "nb", 
    n_layers=2, 
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

## comment if vae model already trained
vae.train()
ref_path = 'ref_model/'
vae.save(ref_path, overwrite=True)

## Uncomment if vae model already trained
#vae =  sca.models.SCVI.load("ref_model/", adata=source_adata)

scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
print("Labelled Indices: ", len(scanvae._labeled_indices))
print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))

# setup and train SCANVI model on reference
print("Train scanvi model")
## comment if scanvae model already trained
scanvae.train(max_epochs=100)
scanvae.save(prefix+"scanvi_model_ref/", overwrite=True) #, save_anndata=True)

## Uncomment if scanvae model already trained
#scanvae =  sca.models.SCANVI.load(prefix+"scanvi_model_ref/", adata=source_adata)

# Setup query adata
target_adata.obs['orig_cell_types'] = target_adata.obs['celltypes'].copy()
target_adata.obs[cell_type_key] = scanvae.unlabeled_category_

## train query adata on SCANVI model
## Comment if model already trained
model = sca.models.SCANVI.load_query_data(
    target_adata,
    prefix+"/scanvi_model_ref/",
    freeze_dropout = True,
)
model._unlabeled_indices = np.arange(target_adata.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))

## Comment if model already trained
model.train(
    max_epochs=200,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
    n_samples_per_label=100
)

# Uncomment if model already trained
#model =  sca.models.SCANVI.load("surgery_model/", adata=source_adata)

predicted = model.predict(adata=target_adata)
predicted_df = pd.DataFrame(predicted, index = target_adata.obs.index)
predicted_df.to_csv(prefix+"/predicted_celltypes.csv")

surgery_path = '/surgery_model'
model.save(prefix+surgery_path, overwrite=True)

# get latent space of reference + query
print("Save latent_representation.csv")
adata_full = source_adata.concatenate(target_adata)
latent_representation = model.get_latent_representation(adata=adata_full)
TOTALVI_LATENT_KEY = "X_totalVI"
adata_full.obsm[TOTALVI_LATENT_KEY] = latent_representation
latent_df = pd.DataFrame(latent_representation, index = adata_full.obs.index)
latent_df.to_csv(prefix + surgery_path+"/latent.csv", index=True)

print("Save umap.csv")
sc.pp.neighbors(adata_full, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(adata_full, min_dist=0.4)
umap_df = pd.DataFrame(adata_full.obsm['X_umap'], index = adata_full.obs.index)
umap_df.to_csv(prefix + surgery_path+"/umap_python.csv", index=True)

print("Save adata_full metadata.csv")
adata_full_metadata = adata_full.obs
adata_full_metadata.to_csv(prefix+"/adata_full_metadata.csv", index=True)


