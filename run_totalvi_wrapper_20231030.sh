#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=gpu 
#SBATCH --gres=gpu:1 
#SBATCH --mem=180GB #184GB
#SBATCH -J totalvi              
#SBATCH -o totalvi_igt1_56_20231030_allgenes.log
#SBATCH -t 24:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch run_totalvi_wrapper_20231030.sh

cd /project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git

module load python/anaconda-2022.05 
module load hdf5/1.12.0 #for BPCells but also others

source activate /project/jfkfloor2/zemmourlab/david/envs/scvi

#Rscript run_totalvi.R /scratch/midway3/zemmour/integration immgent_20231030.Rds totalvi_igt1_56_20231030_allgenes

#python run_totalvi_denoising.py /project/jfkfloor2/zemmourlab/david/immgent/analysis/integration totalvi_igt1_48_20230706_allgenes 

Rscript run_totalvi_plots.R /scratch/midway3/zemmour/integration totalvi_igt1_56_20231030_allgenes


