#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3 
#SBATCH --gres=gpu:1 
#SBATCH --mem=96GB #184GB
#SBATCH -J expimap              
#SBATCH -o expimap_20240703_4.log
#SBATCH -t 12:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch run_totalvi_wrapper_20231030.sh

module load python/anaconda-2022.05 

source activate /project/zemmour/david/envs/scarches_test9

cd /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
working_dir=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg
path_to_anndata=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg/adata_forexpimap.h5ad
prefix=expimap/20240703_4
batchkey=IGT
path_to_signatures=/project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT1_56/Treg/m2.cp.reactome.v2023.2.Mm.symbols.gmt
hvg=500
alpha_kl=0.1
alpha=0.1

cd $working_dir

python $SCRIPT_DIR/run_expimap.py --working_dir=$working_dir --path_to_anndata=$path_to_anndata --prefix=$prefix --batchkey=$batchkey --path_to_signatures=$path_to_signatures --hvg=$hvg --alpha_kl=$alpha_kl --alpha=$alpha
