#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=caslake 
#SBATCH --nodes=1 
#SBATCH --mem=96GB
#SBATCH -J mde_plots_%j            
#SBATCH -o mde_plots_%j.log
#SBATCH -t 12:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/mde_plots_wrapper.sh

module load python/anaconda-2022.05 
source $(conda info --base)/etc/profile.d/conda.sh
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13
conda activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
path_to_wd=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/

cd $path_to_wd

# path_to_seurat_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/igt1_96_withtotalvi20250513_clean.Rds
# output_dir=MDE_InEachSample
# Rscript $SCRIPT_DIR/mde_plots_per_sample.R $path_to_seurat_object $output_dir

path_to_seurat_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/igt1_96_withtotalvi20250710_clean.Rds
output_dir=MDE_InEachSample_20250712
Rscript $SCRIPT_DIR/mde_plots_per_sample.R $path_to_seurat_object $output_dir

