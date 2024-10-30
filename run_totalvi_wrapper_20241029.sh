#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=128GB #184GB
#SBATCH -J totalvi              
#SBATCH -o totalvi_20241029_fitleredGenesProt.log
#SBATCH -t 24:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/run_totalvi_wrapper_20241029.sh

module load python

source activate /project/zemmour/david/envs/scvi120_20241008

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
working_dir=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/totalvi_20241006
path_to_mudata=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/igt1_96_20241006_noTCRRiboMitoLncgenes_TCRCongenicproteins.h5mu
prefix=totalvi_20241029_fitleredGenesProt
batchkey=IGT
categorical_covariate_keys=''
corrected_counts=False
denoised_data=False

cd $working_dir

python $SCRIPT_DIR/run_totalvi_v2.py --working_dir=$working_dir --path_to_mudata=$path_to_mudata --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data
