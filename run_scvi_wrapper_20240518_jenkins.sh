#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=64GB #184GB
#SBATCH -J Jenkins              
#SBATCH -o Jenkins.log
#SBATCH -t 24:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch wrapper_run_scvi_jenkins.py

module load python/anaconda-2022.05 
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13

source activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13 #load again or error

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
WD=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/integration_forJenkins/

cd $WD

python $SCRIPT_DIR/run_scvi.py  --working_dir=$WD --path_to_adata=Cd4Cd8_adata_norep.h5ad --prefix=RmDataset --batchkey=dataset --latent_key=_X_SCVI

python $SCRIPT_DIR/run_scvi.py  --working_dir=$WD --path_to_adata=Cd4Cd8_adata_norep.h5ad --prefix=RmDatasetTimepoint --batchkey=dataset --confoundings='timepoint' --latent_key=_X_SCVI

python $SCRIPT_DIR/run_scvi.py  --working_dir=$WD --path_to_adata=Cd4Cd8_adata_norep.h5ad --prefix=RmDatasetCelltype --batchkey=dataset --confoundings='cell_type' --latent_key=_X_SCVI

python $SCRIPT_DIR/run_scvi.py  --working_dir=$WD --path_to_adata=Cd4Cd8_adata_norep.h5ad --prefix=RmDatasetTimepointCelltype --batchkey=dataset --confoundings='timepoint,cell_type' --latent_key=_X_SCVI
