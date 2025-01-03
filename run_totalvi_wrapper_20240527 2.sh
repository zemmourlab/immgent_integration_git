#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=128GB #184GB
#SBATCH -J totalvi              
#SBATCH -o totalviv2_20240527.log
#SBATCH -t 24:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch run_totalvi_wrapper_20231030.sh

module load python/anaconda-2022.05 
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13

source activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13 #load again or error

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
working_dir=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/
path_to_mudata=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/export_data/totalvi_igt1_56_20231030_allgenes_mdata.h5mu
prefix=totalvi_igt1_56_allgenes_20240526_igtsampleregressedout
batchkey=IGT
categorical_covariate_keys=sample_id
corrected_counts=True
denoised_data=False

cd $working_dir

python $SCRIPT_DIR/run_totalvi_v2.py --working_dir=$working_dir --path_to_mudata=$path_to_mudata --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data
