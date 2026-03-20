#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=beagle3-bigmem 
#SBATCH --nodes=1 
#SBATCH --mem=96GB
#SBATCH -J limma_%j            
#SBATCH -o limma_%j.log
#SBATCH -t 36:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/limma_wrapper_20260318_CD4_tips.sh

module load python/anaconda-2022.05 
source $(conda info --base)/etc/profile.d/conda.sh
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13
conda activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
path_to_wd=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4

cd $path_to_wd

echo "Sampling seurat object"
path_to_seurat_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4_igt1_96_withtotalvi20260206_clean.Rds
output_dir=DGE_limma/20260318_tips
so_file_name=CD4_igt1_96_withtotalvi20260206_clean_sampled.Rds
# Rscript $SCRIPT_DIR/limma_sample_seuratobject_20260223.R $path_to_seurat_object $output_dir $so_file_name

echo "Making tmm file"
path_to_seurat_object=$output_dir/$so_file_name
tmm_file_name=igt1_96_CD4_20260206_tmm.Rds
# Rscript $SCRIPT_DIR/limma_make_tmm_20241217.R $path_to_seurat_object $output_dir $tmm_file_name

echo "Fitting"
path_to_tmm_object=$output_dir/$tmm_file_name
fit_file_name=igt1_96_CD4_20260206_fit.Rds
# Rscript $SCRIPT_DIR/limma_fit_20260223_level2.IGTHT.R $path_to_seurat_object $path_to_tmm_object $output_dir $fit_file_name

echo "Contrasts in Resting"
# path_to_fit_object=$output_dir/$fit_file_name
path_to_tmm_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/DGE_limma/20260206/igt1_96_CD4_20260206_tmm.Rds
path_to_fit_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/DGE_limma/20260206/igt1_96_CD4_20260206_fit.Rds
prefix_file_name=tips
Rscript $SCRIPT_DIR/limma_contrasts_20260318_CD4_tips.R $path_to_seurat_object $path_to_tmm_object $path_to_fit_object $output_dir $prefix_file_name

