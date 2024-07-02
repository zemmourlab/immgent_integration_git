#!/project/zemmour/david/envs/scvi_20240315/bin/python
"""This script run TOTALVI from a Mudata object specifying batch_key and categorical_covariate_keys if needed, """
#author: David Zemmour
#date: 05/27/2024
#run_expimap.py [cwd] [path to mudata .h5mu] [prefix] [batchkey] [] [corrected_counts] [denoised_data]

import warnings; warnings.simplefilter('ignore')
import argparse
parser = argparse.ArgumentParser(description="run ExpiMap from a MuData object specifying ")
parser.add_argument('--working_dir', help='Working directory')
parser.add_argument('--path_to_mudata', help='Path to MuData object')
parser.add_argument('--path_to_signatures', help='Path to GMT signature file')
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
path_to_signatures = args.path_to_signatures
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

# working_dir = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg'
# path_to_mudata = '/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg/totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu'
# prefix = 'expimap_20240628'
# batchkey = 'IGT'
# path_to_signatures = /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg/m2.cp.reactome.v2023.2.Mm.symbols.gmt'


print(f"Working Directory: {working_dir}")
print(f"Path to AnnData: {path_to_mudata}")
print(f"Path to Signture gmt file: {path_to_signatures}")
print(f"Prefix: {prefix}")
print(f"Batch Key: {batchkey}")
print(f"categorical_covariate_keys: {categorical_covariate_keys}")
print(f"corrected_counts: {corrected_counts}")
print(f"denoised_data: {denoised_data}")
#print(f"Latent Key: {latent_key}")

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

os.makedirs(prefix)
os.chdir(working_dir)

print("Reading MuData")
mdata = mu.read(path_to_mudata)

print("Converting to AnnData")
adata = mdata.mod['RNA']
adata.obs['sample_id']

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

#sc.pp.highly_variable_genes(adata, n_top_genes=500, batch_key = 'IGT', subset=True) skipping for now

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

intr_cvae.save(prefix)

print('Plot the latent space of the reference.')

MEAN = False

adata.obsm['X_cvae'] = intr_cvae.get_latent(mean=MEAN, only_active=True)

sc.pp.neighbors(adata, use_rep='X_cvae')
sc.tl.umap(adata)


######
print('Plot the latent space of the reference.')
adata.write(prefix+'/intr_cvae/adata.h5ad')

sc.pl.umap(adata, color=['IGT', 'organ_simplified',], frameon=False, legend_fontsize = )

fig = plt.gcf()  # Get current figure
ax = plt.gca()   # Get current axis

# Adjust the legend properties
legend = ax.get_legend()
if legend:
    legend.set_title_fontsize(12)  # Example: set title font size
    legend.set_fontsize(1)        # Example: set legend font size

# Show or save the modified plot
plt.show()

```

```{python}

adata = mu.read(prefix+'/intr_cvae/adata.h5ad') #totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu") #mu.read("totalvi_igt1_56_allgenes_Treg_20240306/adata.h5mu
intr_cvae = sca.load
intr_cvae = sca.models.EXPIMAP.load(prefix+'/intr_cvae')

adata.obsm['X_cvae'].shape
adata.varm['I'].shape

intr_cvae.latent_directions(adata = adata)

adata.uns['directions']

#Do gene set enrichment test for condition in reference + query using Bayes Factors.

intr_cvae.latent_enrich(groups='organ_simplified', comparison='rest', use_directions=True, adata = adata)
fig = sca.plotting.plot_abs_bfs(adata, yt_step=0.8, scale_y=2.5, fontsize=7)
fig.set_size_inches(20,20)
fig.savefig(prefix+'/intr_cvae/reactome_enrichement_organsimplified_vs_rest.pdf')




```





