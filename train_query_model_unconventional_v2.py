#!/usr/bin/env python3


import pandas as pd
import mudata as mu
import anndata as AnnData
from scipy.sparse import csc_matrix
import scvi
import numpy as np
import scanpy
import matplotlib.pyplot as plt
import mudata
from scipy.io import mmread
from itertools import chain
import torch

# read in created matrix
matrix_query = mmread("matrix/matrix.mtx")
genes = pd.read_csv("matrix/genes.tsv", delimiter = "\t",header=None)
cells_query = pd.read_csv("matrix/barcodes.tsv",header=None)
proteins = pd.read_csv("../proteins.tsv",header=None)
genes = genes.drop(columns=0)
genes = genes.values.astype(str)
cells_query = cells_query.values.astype(str)

csc_matrix_query = csc_matrix(matrix_query)
csc_matrix_query = np.transpose(csc_matrix_query)
adata_q = AnnData.AnnData(csc_matrix_query)


adata_q.var_names = genes.astype(str)
adata_q.obs_names = cells_query.astype(str)

# make protein zeros matrix
protein_data = np.zeros((len(cells_query), 180))
protein_data = protein_data.astype(int)
#protein_data = csc_matrix(protein_data)
#protein_data.index = proteins.astype(str)
#protein_data.columns = cells_query.astype(str)


# read in mudata
adata = mu.read("../model/adata.h5mu") #load this one (results of this block)
print(adata)
adata.mod["RNA"].layers["counts"] = adata.mod["RNA"].X.copy()

# add annotation level 1 and 2, then subset
annotation_level1 = pd.read_csv("../annotation_level1.csv", delimiter = ",", header = None)
annotation_level2 = pd.read_csv("../annotation_level2.csv", delimiter = ",")
annotation_level1 = annotation_level1.drop(columns = 0)
annotation_level2 = annotation_level2.drop(columns = "Unnamed: 0")
adata.obs['annotation_level1'] = annotation_level1.values.astype(str)
adata.obs['annotation_level2'] = annotation_level2.values.astype(str)
adata = adata[adata.obs['annotation_level1'] == "unconventional"]
adata.mod["RNA"].layers["counts"] = adata.mod["RNA"].X.copy()

# setup protein adata
protein_data_adata = adata.mod['protein'].copy()
protein_data_adata.X = np.zeros((len(protein_data_adata), 180))
#protein_data_adata = protein_data_adata[:len(cells_query)]
#protien_adata = protein_data_adata[:len(cells_query), :]
protein_adata = scanpy.AnnData(X=np.array(protein_data_adata[:len(cells_query), :].X))
protein_adata.obs_names = cells_query.astype(str)
protein_adata.var_names = proteins.values.astype(str)

# setup mudata
mdata = mu.MuData({"RNA":adata_q, "protein":protein_adata})
mdata.mod['RNA'].obs_names = cells_query.astype(str)
mdata.mod['RNA'].var_names = genes.astype(str)
mdata.mod['RNA'].layers["counts"] = mdata.mod['RNA'].X.copy()
mdata.mod['protein'].obs_names = cells_query.astype(str)
mdata.mod['protein'].var_names = proteins.values.astype(str)
mdata.obs['batch'] = "query"
mdata.mod['RNA'].obs['batch'] = "query"
mdata.mod['protein'].obs['batch'] = "query"
print(mdata)


# read in model
adata = adata.copy()
model = scvi.model.TOTALVI.load("../model", adata=adata)
model = scvi.model.TOTALVI.load("../model", adata=adata)
model.adata.mod['RNA'].var_names
model.adata.mod['RNA'].obs_names


# setup mudatas
# reference
scvi.model.TOTALVI.setup_mudata(
    adata, 
    rna_layer="counts",
    protein_layer=None,
    batch_key="IGT",
    modalities={
        "rna_layer": "RNA",
        "protein_layer": "protein",
        "batch_key": "RNA",
        "categorical_covariate_keys":"RNA"
    },
)


# Query
scvi.model.TOTALVI.prepare_query_mudata(
  mdata,
  model,
  return_reference_var_names=True
)

#mu.pp.intersect_obs(mudata)
mdata.obs['IGT'] = "query"
mdata.mod['RNA'].obs['IGT'] = "query"
mdata.mod['protein'].obs['IGT'] = "query"

# load query to model
#totalvi_query = scvi.model.TOTALVI.load_query_data(
#    mdata,
#    model
#)

# train query on model
# uncomment if query model not yet trained
totalvi_query.train(200, plan_kwargs={"weight_decay": 0.0})
#torch.save(totalvi_query,"N:/CBDM_Lab/Odhran/scVerse/ImmgenT_Workshop/CD4/BHLHE40/model.pt")
totalvi_query.save("model/model.pt")

# uncomment if model already trained
#totalvi_query = scvi.model.TOTALVI.load("model/model.pt/", adata=mdata)
#totalvi_query = scvi.model.TOTALVI.load("model/model.pt/", adata=mdata)


# get latent representation etc
# plot and visualise integration
# concat just RNA adatas? - for plotting latent space
TOTALVI_LATENT_KEY = "X_totalVI"
mdata.obsm[TOTALVI_LATENT_KEY] = totalvi_query.get_latent_representation()

scanpy.pp.neighbors(mdata, use_rep=TOTALVI_LATENT_KEY)
scanpy.tl.umap(mdata, min_dist=0.4)

query_metadata = pd.read_csv("cell_metadata.csv", delimiter = ";")
query_metadata.index = mdata.obs_names.astype(str)
mdata.obs['annotation_level2'] = query_metadata["celltypes"].values
mdata.obs['organ'] = query_metadata["organ"].values

scanpy.pl.umap(
    mdata,
    color=["annotation_level2"],
    frameon=False,
    ncols=1,
    title="Query",
    save = "Query_test_celltypes.png"
)



adata_query_RNA = mdata.mod['RNA'].copy()
adata_query_RNA.write_h5ad("query_RNA_adata.h5ad")
adata_query_protein = mdata.mod['protein'].copy()
adata_query_protein.write_h5ad("query_protein_adata.h5ad")
latent_query = pd.DataFrame(mdata.obsm["X_totalVI"])
latent_query.to_csv("latent_query.csv")

adata_query = mdata.mod['RNA'].copy()
adata_query.obsm[TOTALVI_LATENT_KEY] = mdata.obsm[TOTALVI_LATENT_KEY]

#adata_ref = adata.mod['RNA'].copy()
#adata_ref.obsm[TOTALVI_LATENT_KEY] = latent_ref

#full = AnnData.concat([adata_ref, adata_query], label = 'batch', keys=['reference','query'])

#latent_ref = pd.read_csv("../latent.csv")
#latent_ref.index = adata.mod['RNA'].obs_names.astype(str)
#latent_ref['annotation_level1'] = annotation_level1.values.astype(str)
#latent_ref = latent_ref[latent_ref["annotation_level1"] == "unconventional"]
#latent_ref = latent_ref.drop(columns = "Unnamed: 0")
#latent_ref = latent_ref.drop(columns = "annotation_level1")
#latent_ref.index = adata.obs_names.astype(str)
#adata.obsm[TOTALVI_LATENT_KEY] = latent_ref

adata.write_h5mu("ref_mudata.h5mu")
#mdata.write_h5mu("query_mudata.h5mu")

mdata.mod['RNA'].var_names = adata.mod['RNA'].var_names
mdata.mod['protein'].var_names = adata.mod['protein'].var_names
mdata.mod['RNA'].obs_names = cells_query.tolist()
mdata.mod['protein'].obs_names = cells_query.tolist()
mdata.obs['batch'] = "query"
adata.obs['batch'] = "reference"
mdata.obs['annotation_level2'] = query_metadata['celltypes'].values
mdata.mod['RNA'].obs['annotation_level2'] = query_metadata['celltypes'].values

full = mu.concat([adata,mdata], label = 'batch', keys=['reference','query'])

full.obsm[TOTALVI_LATENT_KEY] = totalvi_query.get_latent_representation(full)
full.obs['batch']
full.obs['batch'][len(adata.obs_names):] = "query"
full.obs['annotation_level2'][len(adata.obs_names):] = query_metadata['celltypes'].values
#adata.obs['annotation_level2'].values + query_metadata['celltypes'].values
full.obs['clusters/celltypes'] = "NA"
clusters_celltypes = pd.concat([adata.obs['annotation_level2'],query_metadata['celltypes']])
full.obs['clusters/celltypes'] = clusters_celltypes.astype(str)


scanpy.pp.neighbors(full, use_rep=TOTALVI_LATENT_KEY)
scanpy.tl.umap(full, min_dist=0.4)
scanpy.pl.umap(
    full,
    color=["batch"],
    frameon=False,
    ncols=1,
    title="Query & reference",
    save = "Query_Ref_test_full.png"
)

scanpy.pl.umap(
    full,
    color=["clusters/celltypes"],
    frameon=False,
    ncols=1,
    title="Query & reference",
    save = "Query_Ref_test_celltypes_full.png"
)

full.write_h5mu("full_mudata.h5mu")

