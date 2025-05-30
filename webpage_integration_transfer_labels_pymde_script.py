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
#mdata = mu.read(path_to_mudata) #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
#print(mdata)
#mdata = AnnData.read(working_dir+"/MERGED_RNA_Manual_GSE182509_samples_processed.h5ad") #load this one (results of this block)
ImmgenT_mdata = AnnData.read(path_to_ImmgenT)
ImmgenT_mdata.X = ImmgenT_mdata.X.copy()
ImmgenT_mdata.layers["counts"] = ImmgenT_mdata.X.copy()
print(ImmgenT_mdata)
query_mdata = AnnData.read(path_to_query)
query_mdata.X = query_mdata.X.copy()
query_mdata.layers["counts"] = query_mdata.X.copy()
print(query_mdata)

#full.mod['RNA'].obs_names = full_RNA.obs_names.astype(str)
#full.mod['RNA'].var_names = genes['x'].astype(str)
#full.mod['ADT'].obs_names = full_ADT.obs_names.astype(str)
#full.mod['ADT'].var_names = proteins.astype(str)

#ImmgenT_IGTHT = pd.read_csv(working_dir+"ImmgenT_IGTHT.csv", delimiter = ",")
#ImmgenT_IGTHT.index = ImmgenT_IGTHT['Unnamed: 0'].values
#ImmgenT_IGTHT = ImmgenT_IGTHT.drop(columns = "Unnamed: 0")
#ImmgenT_mdata.obs['IGTHT'] = ImmgenT_IGTHT["IGTHT"].astype(str)
#ImmgenT_mdata.obs['IGT'] = ImmgenT_IGTHT["IGT"].astype(str)


#query_IGTHT = pd.read_csv(working_dir+"query_IGTHT.csv", delimiter = ",")
#query_IGTHT.index = query_IGTHT['Unnamed: 0'].values
#query_IGTHT = query_IGTHT.drop(columns = "Unnamed: 0")
#query_IGTHT['Unknown'] = "Unknown"

print("Concat anndata - keep only common genes")
mdata = AnnData.concat(
    [ImmgenT_mdata, query_mdata],
    axis=0,
    join="inner",
    label="origin",
    keys=["ImmgenT", "query"]
)

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
mdata.obs['level2.group'] = "Unknown"
mdata.obs.loc[annotation_level2.index, 'level1'] = annotation_level2['level1'].values
mdata.obs.loc[annotation_level2.index, 'level2'] = annotation_level2['level2'].values
mdata.obs.loc[annotation_level2.index, 'level2.group'] = annotation_level2['level2.group'].values


# Comment if already trained
scvi.model.SCVI.setup_anndata(mdata, layer = "counts", batch_key = batchkey, categorical_covariate_keys = categorical_covariate_keys)
scvi_model = scvi.model.SCVI(mdata, n_latent=30, n_layers=2)
scvi_model.train()
scvi_model.save(prefix+"/scvi_model/") #, save_anndata=True)

SCANVI_LATENT_KEY = "X_scANVI"
print("Save latent_representation.csv")
latent_representation = scvi_model.get_latent_representation()
TOTALVI_LATENT_KEY = "X_totalVI"
##mdata.mod['RNA'].obsm[TOTALVI_LATENT_KEY] = latent_representation
mdata.obsm[TOTALVI_LATENT_KEY] = latent_representation
##latent_df = pd.DataFrame(latent_representation, index = mdata.mod['RNA'].obs.index)
latent_df = pd.DataFrame(latent_representation, index = mdata.obs.index)

latent_df.to_csv(prefix+"/latent.csv", index=True)

print("Save umap.csv")
sc.pp.neighbors(mdata, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(mdata, min_dist=0.4)
##umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.mod['RNA'].obs.index)
umap_df = pd.DataFrame(mdata.obsm['X_umap'], index = mdata.obs.index)
umap_df.to_csv(prefix+"/umap_python.csv", index=True)

print("Save Anndata")
mdata.write_h5ad(prefix+"/adata_RNA.h5ad")

#mdata = AnnData.read(prefix+"_orig/adata_RNA.h5ad")
#latent_df = pd.read_csv(prefix+"_orig/latent.csv", index_col=0)

print("Transfer labels using Annoy - level1 and 2")
# Define latent embedding and masks
latent = latent_df
##mdata.obs["is_ref"] = mdata.obs["origin"]
ref_mask = mdata.obs["origin"] == "ImmgenT"
query_mask = mdata.obs["origin"] == "query"
ref_embeddings = latent[ref_mask.values]
query_embeddings = latent[query_mask.values]

# Build Annoy index on ref
n_neighbors = 30
n_dims = ref_embeddings.shape[1]
annoy_index = AnnoyIndex(n_dims, "angular")
annoy_index.set_seed(0) # ensure reproducibility

##for i, vec in enumerate(ref_embeddings):
##    annoy_index.add_item(i, vec)
##annoy_index.build(10)

for i in range(ref_embeddings.shape[0]):
    annoy_index.add_item(i, ref_embeddings.iloc[i].tolist())

annoy_index.build(10)


# Find neighbors for each query cell
##nn_indices = np.vstack([
##    annoy_index.get_nns_by_vector(vec, n_neighbors) for vec in query_embeddings
##])

##nn_indices = np.vstack([
##    annoy_index.get_nns_by_vector(vec.tolist(), n_neighbors)
##    for vec in query_embeddings.to_numpy()
##])

neighbors_list = []
for vec in query_embeddings.to_numpy():
    neighbors = annoy_index.get_nns_by_vector(vec.tolist(), n_neighbors)
    neighbors_list.append(neighbors)

nn_indices = np.vstack(neighbors_list)
print(nn_indices)
print(nn_indices.shape)

## Loop over annotation levels
for level in ["level1", "level2", "level2.group"]:
    ## Extract true ref labels
    label_column = f"{level}"
    new_column = f"{level}_transfer_labels"
    
    ## Initialize with original annotations
    mdata.obs[new_column] = mdata.obs[label_column].astype(str)

    ## Transfer labels via majority vote
    labels_ref = mdata.obs.loc[ref_mask, label_column].astype(str).values
    ##labels_ref = mdata.obs.loc[mdata.obs['origin'] == "ImmgenT", label_column].astype(str).values
    ##labels_ref = mdata[ref_mask].obs[label_column].astype(str).values
    
    transfer_labels = []
    for neighbors in nn_indices:
        neighbor_labels = labels_ref[neighbors]
        most_common = Counter(neighbor_labels).most_common(1)[0][0]
        transfer_labels.append(most_common)
    
    ##transfer_labels = [
    ##    Counter(labels_ref[neighbors]).most_common(1)[0][0]
    ##    for neighbors in nn_indices
    ##]
    ##for neighbors in nn_indices:
    ##    neighbor_labels = labels_ref[neighbors]  # These indices are relative to ref_embeddings
    ##    most_common = Counter(neighbor_labels).most_common(1)[0][0]
    ##    transfer_labels.append(most_common)
    
    ## Assign transferred labels to query cells
    mdata.obs.loc[query_mask, new_column] = transfer_labels
    
    ## Get the indices of the query cells
    ##query_indices = mdata.obs.index[query_mask]
    ##mdata.obs.loc[query_indices, new_column] = pd.Series(transfer_labels, index=query_indices)
    ## Create a pandas Series with proper index alignment
    ##transfer_series = pd.Series(transfer_labels, index=query_indices)

    ## Assign transferred labels to the correct cells
    ##mdata.obs.loc[query_indices, new_column] = transfer_series

output_file = mdata.obs[['level1_transfer_labels','level2_transfer_labels','level2.group_transfer_labels']]
output_file.index = mdata.obs.index.copy()
output_file.to_csv(prefix+"/output_annotations.csv", index=True)

print("run pyMDE - first on all data, second on subgroups")
print("pyMDE on all data")
# Load the full MDE reference embedding
mde_ref_embedding = pd.read_csv(mde_ref_file, index_col=0)
mde_ref_embedding = mde_ref_embedding.loc[mde_ref_embedding.index.intersection(latent_df.index)]

# level1 mde
rownames = mde_ref_embedding.index.tolist()  # Convert index to list
index_positions = [latent_df.index.get_loc(item) for item in rownames if item in latent_df.index]
print(index_positions)
print(len(index_positions))
print(len(mde_ref_embedding.index))
rownames_tensor = torch.tensor(index_positions, device='cpu') # , device='cuda:0')

print("pymde.Anchored")
anchor_constraint = pymde.Anchored(
    anchors=rownames_tensor,
    values=torch.tensor(mde_ref_embedding.values, dtype=torch.float32, device='cpu'), # device='cuda:0),' #device='cpu' 'cuda:0'
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
    group_cells = mdata.obs_names[mdata.obs["level1_transfer_labels"] == annotation]
    
    # Subset the totalVI embedding to group
    group_integrated = latent_df.loc[latent_df.index.intersection(group_cells)]
    
    # Also subset reference anchors to those in this group
    group_ref_embedding = mde_ref_embedding.loc[mde_ref_embedding.index.intersection(group_integrated.index)]
    
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
user_output_file = output_file.loc[query_mdata.obs.index, :] 
output_file.to_csv(prefix+"/output_file.csv", index=True)
user_output_file.index = query_mdata.obs.index
user_output_file.to_csv(prefix+"/user_output_file.csv", index=True)

