#!/bin/bash
#SBATCH -A pi-zemmour 
#SBATCH --partition=gpu 
#SBATCH --gres=gpu:1 
#SBATCH --mem=64GB
#SBATCH -J totalvi              
#SBATCH -o totalvi_igt12_56_20231119_allgenes.log
#SBATCH -t 12:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch run_totalvi_wrapper_20231030.sh

cd /project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git

module load python/anaconda-2022.05 
module load hdf5/1.12.0 #for BPCells but also others

source activate /project/jfkfloor2/zemmourlab/david/envs/scvi

Rscript run_totalvi.R /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT22_56_shanelle immgent_IGT22_56_shanelle.Rds totalvi_igt22_56_shanelle_20231119

#python run_totalvi_denoising.py /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration totalvi_igt1_48_20230706_allgenes 

Rscript run_totalvi_plots.R /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration/IGT22_56_shanelle totalvi_igt22_56_shanelle_20231119


