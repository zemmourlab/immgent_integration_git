#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=64GB #184GB
#SBATCH -J totalvi              
#SBATCH -o totalvi_igt1_56_allgenes_rmorgan2.log
#SBATCH -t 48:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch run_totalvi_wrapper.sh

module load python/anaconda-2022.05 
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13

source activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13 #load again or error

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
CWD=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56
cd $CWD

# python $SCRIPT_DIR/run_totalvi.py $CWD $CWD/export_data/totalvi_igt1_56_20231030_allgenes_mdata.h5mu totalvi_igt1_56_20231030_allgenes_organregressed IGT organ_simplified X_totalVI_rmorgan

python $SCRIPT_DIR/run_totalvi_denoising_v2.py $CWD totalvi_igt1_56_20231030_allgenes_organregressed
