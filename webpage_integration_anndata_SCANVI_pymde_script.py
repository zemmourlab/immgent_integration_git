#!/usr/bin/env python3
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 10/08/2024
#run_totalvi_v2.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [categorical_covariate_keys] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run TOTALVI from a MuData object specifying batch_key and categorical_covariate_keys if needed, can also return corrected counts and denoised data")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_ImmgenT', help='Path to ImmgenT MuData object')
parser.add_argument('--path_to_query', help='Path to query MuData object')
parser.add_argument('--path_to_anndata', help='Path to combined AnnData object')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')
parser.add_argument('--categorical_covariate_keys', default=None, help='categorical_covariate_keys variables (default: None)')
parser.add_argument('--corrected_counts', default=False, help='Returns corrected counts, aka posterior_predictive_sample() (default: False)')
parser.add_argument('--denoised_data', default=False, help='Returns denoised data, aka get_normalized_expression()  (default: False)')
parser.add_argument('--mde_ref_file', help='mde_ref_file  (default: None)')
#parser.add_argument('--latent_key', help='Key for latent space')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_ImmgenT = args.path_to_ImmgenT
path_to_query = args.path_to_query
path_to_anndata = args.path_to_anndata
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
mde_ref_file = args.mde_ref_file
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
print(f"Path to ImmgenT AnnData: {path_to_ImmgenT}")
print(f"Path to query AnnData: {path_to_query}")
print(f"Path to combined AnnData: {path_to_anndata}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"categorical_covariate_keys: {categorical_covariate_keys}")
print(f"corrected_counts: {corrected_counts}")
print(f"denoised_data: {denoised_data}")
print(f"mde_ref_file: {mde_ref_file}")
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
import pymde #to run MDE
from annoy import AnnoyIndex
from collections import Counter

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
##mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
##print(mdata)
mdata = AnnData.read(path_to_anndata) #load this one (results of this block)
mdata.raw = None
#ImmgenT_mdata = AnnData.read(path_to_ImmgenT)
#ImmgenT_mdata.X = ImmgenT_mdata.X.copy()
#ImmgenT_mdata.layers["counts"] = ImmgenT_mdata.X.copy()
#print(ImmgenT_mdata)
## Read in query Anndata object
#query_mdata = AnnData.read(path_to_query)
#query_mdata.X = query_mdata.X.copy()
#query_mdata.layers["counts"] = query_mdata.X.copy()
#print(query_mdata)

##full.mod['RNA'].obs_names = full_RNA.obs_names.astype(str)
##full.mod['RNA'].var_names = genes['x'].astype(str)
##full.mod['ADT'].obs_names = full_ADT.obs_names.astype(str)
##full.mod['ADT'].var_names = proteins.astype(str)

##ImmgenT_IGTHT = pd.read_csv(working_dir+"ImmgenT_IGTHT.csv", delimiter = ",")
##ImmgenT_IGTHT.index = ImmgenT_IGTHT['Unnamed: 0'].values
##ImmgenT_IGTHT = ImmgenT_IGTHT.drop(columns = "Unnamed: 0")
##ImmgenT_mdata.obs['IGTHT'] = ImmgenT_IGTHT["IGTHT"].astype(str)
##ImmgenT_mdata.obs['IGT'] = ImmgenT_IGTHT["IGT"].astype(str)


##query_IGTHT = pd.read_csv(working_dir+"query_IGTHT.csv", delimiter = ",")
##query_IGTHT.index = query_IGTHT['Unnamed: 0'].values
##query_IGTHT = query_IGTHT.drop(columns = "Unnamed: 0")
##query_IGTHT['Unknown'] = "Unknown"

#print("Concat anndata - keep only common genes")
#mdata = AnnData.concat(
#    [ImmgenT_mdata, query_mdata],
#    axis=0,
#    join="inner",
#    label="origin",
#    keys=["ImmgenT", "query"]
#)

print("Creating RNA layer counts")
mdata.layers["counts"] = mdata.X.copy()

##mdata.obs['IGTHT'] = pd.concat([ImmgenT_IGTHT["IGTHT"],query_IGTHT["IGTHT"]]).astype(str)
##mdata.obs['IGT'] = pd.concat([ImmgenT_IGTHT["IGT"],query_IGTHT["IGT"]]).astype(str)
mdata.layers["counts"] = mdata.X.copy()
mdata.obs_names = mdata.obs_names.astype(str)
mdata.var_names = mdata.var_names.astype(str)
mdata.X = mdata.X.copy()
##mdata.var = mdata.var.reset_index(drop=True)
##mdata.mod['RNA'].var = mdata.mod['RNA'].var.reset_index(drop=True)
##mdata.mod['protein'].var = mdata.mod['protein'].var.reset_index(drop=True)
print(mdata)
##mdata.update()

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
        ##print(mdata.mod["RNA"].obs[c].head(10))
        ##mdata.mod["RNA"].obs[c] = mdata.mod["RNA"].obs[c].str.replace('.', '', regex=False)
        ##print(mdata.mod["RNA"].obs[c].head(10))

print("Train new SCVI model")
annotation_level2 = pd.read_csv("/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_annotation_table.csv", delimiter = ",", index_col=0)
annotation_level2 = annotation_level2.loc[annotation_level2.index.intersection(mdata.obs.index)]
mdata.obs['level1'] = "Unknown"
mdata.obs['level2'] = "Unknown"
#mdata.obs['level2.group'] = "Unknown"
mdata.obs.loc[annotation_level2.index, 'level1'] = annotation_level2['level1'].values
mdata.obs.loc[annotation_level2.index, 'level2'] = annotation_level2['level2'].values
#mdata.obs.loc[annotation_level2.index, 'level2.group'] = annotation_level2['level2.group'].values


# Comment if already trained
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
scvi_model = scvi.model.SCVI(mdata, n_latent=30, n_layers=2)
scvi_model.train()
scvi_model.save(prefix+"/scvi_model/") #, save_anndata=True)

#SCANVI_LATENT_KEY = "X_scANVI"
print("Save latent_representation.csv")
latent_representation = scvi_model.get_latent_representation()
SCVI_LATENT_KEY = "X_scVI"
##mdata.mod['RNA'].obsm[SCVI_LATENT_KEY] = latent_representation
mdata.obsm[SCVI_LATENT_KEY] = latent_representation
##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
latent_df = pd.DataFrame(latent_representation, index = mdata.obs.index)

latent_df.to_csv(prefix+"/latent.csv", index=True)

print("Save umap.csv")
sc.pp.neighbors(mdata, use_rep=SCVI_LATENT_KEY)
sc.tl.umap(mdata, min_dist=0.4)
##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.obs.index)
umap_df.to_csv(prefix+"/umap_python.csv", index=True)

print("Save Anndata")
mdata.write_h5ad(prefix+"/adata_RNA.h5ad")

#mdata = AnnData.read(prefix+"_orig/adata_RNA.h5ad")
#latent_df = pd.read_csv(prefix+"_orig/latent.csv", index_col=0)

print("Predict level1 and 2 annotations with SCANVI")
## Define ref and query masks (future use)
ref_mask = mdata.obs["origin"] == "ImmgenT"
query_mask = mdata.obs["origin"] == "query"

## level1
## Training scanvi model on scvi model
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys, labels_key="level1")
level1_model = scvi.model.SCANVI.from_scvi_model(scvi_model, "Unknown")
level1_model.train(25)
level1_model.save(prefix+"/scanvi_level1_model/") #, save_anndata=True)

## level2
## Training scanvi model on scvi model
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys, labels_key="level2")
level2_model = scvi.model.SCANVI.from_scvi_model(scvi_model, "Unknown")
level2_model.train(25)
level2_model.save(prefix+"/scanvi_level2__model/") #, save_anndata=True)

## Predictions and scores - create output file
## level1
#SCVI_LATENT_KEY = "X_SCVI"
LEVEL1_SCANVI_LATENT_KEY = "level1_X_scANVI"
LEVEL1_SCANVI_PREDICTIONS_KEY = "level1_C_scANVI"

#mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(mdata)
mdata.obsm[LEEVL1_SCANVI_LATENT_KEY] = level1_model.get_latent_representation(mdata)
mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY]= level1_model.predict(mdata)
output_file = pd.DataFrame(mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY], index = mdata.obs.index)

# Get posterior probabilities for all labels
level1_probs = model.predict(mdata, soft=True)
# Get max probability per cell (i.e., model confidence)
level1_confidence = probs.max(axis=1)
# Add to AnnData
output_file["level1_scanvi_confidence"] = level1_confidence
#output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

# Add final annotation  with unclear below a threshold
confidence_threshold = 0.95
output_file["level1_final"] = output_file[LEVEL1_SCANVI_PREDICTIONS_KEY]
output_file.loc[output_file["level1_scanvi_confidence"] < confidence_threshold, "level1_final"] = "unclear"

## level2
#SCVI_LATENT_KEY = "X_SCVI"
LEVEL2_SCANVI_LATENT_KEY = "level2_X_scANVI"
LEVEL2_SCANVI_PREDICTIONS_KEY = "level2_C_scANVI"

#mdata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation(mdata)
mdata.obsm[LEEVL2_SCANVI_LATENT_KEY] = level2_model.get_latent_representation(mdata)
mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]= level2_model.predict(mdata)
output_file[LEVEL2_SCANVI_PREDICTIONS_KEY] = mdata.obs[LEVEL2_SCANVI_PREDICTIONS_KEY]

## Get posterior probabilities for all labels
level2_probs = model.predict(mdata, soft=True)
## Get max probability per cell (i.e., model confidence)
level2_confidence = probs.max(axis=1)
## Add to AnnData
output_file["level2_scanvi_confidence"] = level2_confidence
#output_file.to_csv(prefix_SCANVI+"/predicted_celltypes.csv")

## Add final annotation  with unclear below a threshold
output_file["level2_final"] = output_file[LEVEL2_SCANVI_PREDICTIONS_KEY]
output_file.loc[output_file["level2_scanvi_confidence"] < confidence_threshold, "level2_final"] = "unclear"

# Save user annotations to return to the user
##output_file = mdata.obs[['level1_transfer_labels','level2_transfer_labels']] #,'level2.group_transfer_labels']]
##output_file.index = mdata.obs.index.copy()
##output_file.to_csv(prefix+"/output_annotations.csv", index=True)
user_output_file = output_file.loc[query_mask, :]
user_output_file.index = query_mask.obs.index
user_output_file.to_csv(prefix+"/user_output_file.csv", index=True)



print("run pyMDE - first on all data, second on subgroups")
print("pyMDE on all data")
latent_df = mdata.obsm[LEVEL1_SCANVI_KEY]
# Load the full MDE reference embedding
mde_ref_embedding = pd.read_csv(mde_ref_file, index_col=0)
mde_ref_embedding = mde_ref_embedding.loc[mde_ref_embedding.index.intersection(latent_df.index)]

# Make level1 and level2 mde files
# level1
mde_ref_embedding_level1 = mde_ref_embedding[["level1_MDE1","level1_MDE2"]]
mde_ref_embedding_level1.index = mde_ref_embedding.index.copy()

# level2
mde_ref_embedding_level2 = mde_ref_embedding[["level2_MDE1","level2_MDE2"]]
mde_ref_embedding_level2.index = mde_ref_embedding.index.copy()


# level1 mde
rownames = mde_ref_embedding_level1.index.tolist()  # Convert index to list
index_positions = [latent_df.index.get_loc(item) for item in rownames if item in latent_df.index]
print(index_positions)
print(len(index_positions))
print(len(mde_ref_embedding_level1.index))
rownames_tensor = torch.tensor(index_positions, device='cpu') # , device='cuda:0')

print("pymde.Anchored")
anchor_constraint = pymde.Anchored(
    anchors=rownames_tensor,
    values=torch.tensor(mde_ref_embedding_level1.values, dtype=torch.float32, device='cpu'), # device='cuda:0),' #device='cpu' 'cuda:0'
)

print("pymde.preserve_neighbors")
incremental_mde = pymde.preserve_neighbors(
    torch.tensor(latent_df.values, device='cpu'), # dtype=torch.float32, device='cpu'), #device='cuda:0'), #device='cpu'
    embedding_dim=2,
    constraint=anchor_constraint,
    repulsive_fraction=0.7,
    n_neighbors = 15,
    verbose=True,
    device='cpu'  # 'cuda:0'
)

print("incremental_mde.embed")
incremental_mde.embed(eps=1e-6, verbose=True)
level1_mde_df = pd.DataFrame(incremental_mde.X.cpu().numpy(), index = latent_df.index)

output_file['level1_MDE1'] =  level1_mde_df[0]
output_file['level2_MDE1'] =  level1_mde_df[0]
output_file['level1_MDE2'] =  level1_mde_df[1]
output_file['level2_MDE2'] =  level1_mde_df[1]

print("pyMDE on subgroups")
# Subset to level1 annotations
##level1_annotations = mdata.obs["level1_transfer_labels"].unique()
level1_annotations = ["CD4", "CD8", "Treg", "gdT", "CD8aa", "nonconv", "DN", "DP"]

# Directory to save outputs
##os.makedirs(f"{prefix}/mde_by_level1", exist_ok=True)

for annotation in level1_annotations:
    print(f"\n--- Processing group: {annotation} ---")
    
    # Cells in the current group
    group_cells = mdata.obs_names[mdata.obs[LEVEL1_SCANVI_PREDICTIONS_KEY] == annotation]
    
    # Subset the totalVI embedding to group
    group_integrated = latent_df.loc[latent_df.index.intersection(group_cells)]
    
    # Also subset reference anchors to those in this group
    group_ref_embedding = mde_ref_embedding_level2.loc[mde_ref_embedding_level2.index.intersection(group_integrated.index)]
    
    # Check if we have enough anchors to continue
    if group_ref_embedding.shape[0] < 5:
        print(f"Skipping {annotation} (too few anchors: {group_ref_embedding.shape[0]})")
        continue
    
    if group_integrated.shape[0] < 1:
        print(f"Skipping {annotation} (too few {annotation} cells: {group_integrated.shape[0]})")
        continue


    # Get positions of anchor rows in SCVI
    index_positions = [
        group_integrated.index.get_loc(cell)
        for cell in group_ref_embedding.index if cell in group_integrated.index
    ]
    
    if len(index_positions) != group_ref_embedding.shape[0]:
        print("Mismatch between index positions and anchors — skipping.")
        continue

    # Prepare PyTorch tensors
    rownames_tensor = torch.tensor(index_positions, device='cpu')
    anchor_values = torch.tensor(group_ref_embedding.values, dtype=torch.float32, device='cpu')
    full_tensor = torch.tensor(group_integrated.values, dtype=torch.float32, device='cpu')

    # Anchored MDE constraint
    print("  > Building anchor constraint")
    anchor_constraint = pymde.Anchored(anchors=rownames_tensor, values=anchor_values)

    # Run MDE
    print("  > Running pymde.preserve_neighbors")
    mde_model = pymde.preserve_neighbors(
        full_tensor,
        embedding_dim=2,
        constraint=anchor_constraint,
        repulsive_fraction=0.7,
        n_neighbors=15,
        verbose=True,
        device='cpu'
    )

    print("  > Embedding")
    mde_model.embed(eps=1e-6, verbose=True)

    # Add result to output file
    level2_subgroup_mde_df = pd.DataFrame(mde_model.X.cpu().numpy(), index=group_integrated.index)
    output_file.loc[group_integrated.index, 'level2_MDE1'] = level2_subgroup_mde_df[0].values
    output_file.loc[group_integrated.index, 'level2_MDE2'] = level2_subgroup_mde_df[1].values
    ##out_path = f"{prefix}/mde_by_level1/mde_{annotation.replace(' ', '_')}.csv"
    ##output_df.to_csv(out_path)
    ##print(f"  > Saved to {out_path}")


print("Prepare output files - plotting and annotations for user")
#user_output_file = output_file.loc[query_mask.obs.index, :] 
output_file.to_csv(prefix+"/output_file.csv", index=True)
#user_output_file.index = query_mask.obs.index
#user_output_file.to_csv(prefix+"/user_output_file.csv", index=True)

