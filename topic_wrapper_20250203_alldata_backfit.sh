#!/bin/bash
#SBATCH -A pi-zemmour 
#SBATCH --qos=zemmour
#SBATCH --partition=zemmour-hm
#SBATCH --nodes=1 
#SBATCH --mem=0
#SBATCH -J topic_%j            
#SBATCH -o topic_%j.log
#SBATCH -t 7-00:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/topic_wrapper_20250203_alldata_backfit.sh

module load python/anaconda-2022.05 
source $(conda info --base)/etc/profile.d/conda.sh
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13
conda activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
path_to_wd=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/

cd $path_to_wd

path_to_seurat_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/igt1_96_withtotalvi20250109_clean.Rds
output_dir=topic/flashier20250203_alldata_backfit

Rscript $SCRIPT_DIR/topic_flashier_20250203_withbackfit.R $path_to_seurat_object $output_dir

