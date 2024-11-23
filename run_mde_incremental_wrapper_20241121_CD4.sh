#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=128GB #184GB
#SBATCH -J mde              
#SBATCH -o mde_incremental_20241121_CD4.log
#SBATCH -t 12:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/run_mde_incremental_wrapper_20241121_CD4.sh

module load python

source activate /project/zemmour/david/envs/scvi120_20241008

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
working_dir=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4
# path_to_mudata=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/export_data/igt1_96_20241113_CD8AB.h5mu
prefix=totalvi_20241113_CD4_rmIGTsample
mde_ref_file=totalvi_20241014_CD4_rmIGTsample/mde2.csv
totalvi_integrated_file=totalvi_20241113_CD4_rmIGTsample/latent.csv

cd $working_dir
# python $SCRIPT_DIR/run_mde_incremental.py --working_dir=$working_dir --path_to_mudata=$path_to_mudata --prefix=$prefix --mde_ref_file=$mde_ref_file --totalvi_integrated_file=$totalvi_integrated_file
python $SCRIPT_DIR/run_mde_incremental.py --working_dir=$working_dir --prefix=$prefix --mde_ref_file=$mde_ref_file --totalvi_integrated_file=$totalvi_integrated_file

