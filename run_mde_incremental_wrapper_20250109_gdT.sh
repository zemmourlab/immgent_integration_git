#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=64GB #184GB
#SBATCH -J mde              
#SBATCH -o mde_incremental_20250109_gdT.log
#SBATCH -t 12:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/run_mde_incremental_wrapper_20250109_gdT.sh

module load python
source $(conda info --base)/etc/profile.d/conda.sh

source activate /project/zemmour/david/envs/scvi120_20241008

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
working_dir=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/gdT
prefix=20250109_gdT
mde_ref_file=totalvi_20241201_gdT_rmIGTsample/mde2_nooutliers.csv
totalvi_integrated_file=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/totalvi_20241008_rmIGTsample/latent_gdT.csv
cd $working_dir
python $SCRIPT_DIR/run_mde_incremental.py --working_dir=$working_dir --prefix=$prefix --mde_ref_file=$mde_ref_file --totalvi_integrated_file=$totalvi_integrated_file
