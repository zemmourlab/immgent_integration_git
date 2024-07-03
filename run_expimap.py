#!/project/zemmour/david/envs/scarches_test9/bin/python
"""This script run EXPIMAP from a AnnData object"""
#author: David Zemmour
#date: 07/01/2024
#run_expimap.py 

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run ExpiMap from a AnnData object ")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_anndata', help='Path to AnnData object')
parser.add_argument('--path_to_signatures', help='Path to GMT signature file')
parser.add_argument('--prefix', default='myprefix', 
                    help='Prefix for the output files (default: myprefix)')
parser.add_argument('--batchkey', default=None, help='Batch key for analysis')

print("Arguments")
args = parser.parse_args()
working_dir = args.working_dir
path_to_anndata = args.path_to_anndata
path_to_signatures = args.path_to_signatures
prefix = args.prefix
batchkey = args.batchkey


# working_dir = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg'
# path_to_anndata = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg/adata_forexpimap.h5ad'
# prefix = 'expimap_20240701'
# batchkey = 'IGT'
# path_to_signatures = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg/m2.cp.reactome.v2023.2.Mm.symbols.gmt'


print(f"Working Directory: {working_dir}")
print(f"Path to AnnData: {path_to_anndata}")
print(f"Path to Signture gmt file: {path_to_signatures}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")

print("Importing libraries")
import warnings
warnings.simplefilter(action='ignore')
import scanpy as sc
import torch
import scarches as sca
import numpy as np
import gdown
import mudata as mu
import os
import pandas as pd
import matplotlib.pyplot as plt

print("Global configurations")
# pd.set_option('display.max_rows', 1000)
# pd.set_option('display.max_columns', 1000)
# sc.set_figure_params(figsize=(6, 6), frameon=False)
# scvi.settings.seed = 0  # optional: ensures reproducibility
# #sns.set_theme()
# if torch.cuda.is_available():
#     print("CUDA is available")
#     print("Using device:", torch.cuda.get_device_name())
#     torch.set_float32_matmul_precision("high")

os.chdir(working_dir)

print('Create output directory: {prefix}')
if not os.path.exists(prefix):
    os.makedirs(prefix)
    print(f"Created directory: {prefix}")
else:
    print(f"Directory already exists: {prefix}")

print("Reading AnnData")
#mdata = mu.read(path_to_anndata)
adata = sc.read(path_to_anndata)

print("Converting to AnnData")
#adata = mdata.mod['RNA']
#adata.obs['sample_id']

adata.var_names = adata.var_names.astype(str)
adata.obs_names = adata.obs_names.astype(str)
assert adata.var_names.is_unique, "Variable names are not unique!"
assert adata.obs_names.is_unique, "Observation names are not unique!"

adata.X = adata.layers["counts"].copy() #.X should contain raw counts

print('Reading annotations')
sca.utils.add_annotations(adata, path_to_signatures, min_genes=12, clean=True)
mask = adata.varm['I'].sum(axis=1) > 0
print(f"Number of genes to kept from annotation: {np.where(mask)[0].shape[0]}")

print('Normalization, log1p')
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=500, batch_key = 'IGT', subset=True) #skipping for now

print('Filter out all annotations (terms) with less than 12 genes.')
select_terms = adata.varm['I'].sum(0)>12
adata.uns['terms'] = np.array(adata.uns['terms'])[select_terms].tolist()
adata.varm['I'] = adata.varm['I'][:, select_terms]

print('Filter out genes not present in any of the terms after selection of HVGs.')
adata._inplace_subset_var(adata.varm['I'].sum(1)>0)
adata.var
adata.X = adata.layers["counts"].copy()

print('Create expiMap model')
intr_cvae = sca.models.EXPIMAP(
    adata=adata,
    condition_key='IGT',
    hidden_layer_sizes=[256, 256, 256],
    recon_loss='nb'
)

ALPHA = 0.7
early_stopping_kwargs = {
    "early_stopping_metric": "val_unweighted_loss", # val_unweighted_loss
    "threshold": 0,
    "patience": 50,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}
intr_cvae.train(
    n_epochs=400,
    alpha_epoch_anneal=100,
    alpha=ALPHA,
    alpha_kl=0.5,
    weight_decay=0.,
    early_stopping_kwargs=early_stopping_kwargs,
    use_early_stopping=True,
    monitor_only_val=False,
    seed=2020,
)

intr_cvae.save(prefix+'/model)

print("Save latent data and umap in csv file")
MEAN = False
latent_representation = intr_cvae.get_latent(mean=MEAN, only_active=True)
adata.obsm['X_expimap'] = latent_representation
latent_df = pd.DataFrame(latent_representation, index = adata.obs.index)
latent_df.to_csv(prefix+"/latent.csv", index=True)

sc.pp.neighbors(adata, use_rep='X_expimap')
sc.tl.umap(adata, n_components = 2)
#adata.obsm['X_umap_expimap'] = adata.obsm['X_umap']
#del adata.obsm['X_umap']
umap_df = pd.DataFrame(adata.obsm['X_umap'], index = adata.obs.index, columns = ['umap_1', 'umap_2'])
umap_df.to_csv(prefix+"/umap.csv", index=True)

print('Plot the UMAP')
sc.pl.umap(adata, color=['IGT', 'organ_simplified',], frameon=False, legend_fontsize = 5, show = False)
plt.savefig(prefix+'/umap.pdf')

print('Done')



