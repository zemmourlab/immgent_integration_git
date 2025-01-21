#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=caslake 
#SBATCH --nodes=1 
#SBATCH --mem=96GB
#SBATCH -J limma_%j            
#SBATCH -o limma_%j.log
#SBATCH -t 36:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/limma_wrapper_20250109_gdT.sh

module load python/anaconda-2022.05 
source $(conda info --base)/etc/profile.d/conda.sh
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13
conda activate /project/zemmour/david/envs/scvi_20240315
module load openblas/0.3.13

SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
path_to_wd=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/gdT

cd $path_to_wd

echo "Sampling seurat object"
path_to_seurat_object=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/gdT/igt1_96_gdT_20250109.Rds
output_dir=DGE_limma/20250109
so_file_name=igt1_96_gdT_20250109_sampled_so.Rds
Rscript $SCRIPT_DIR/limma_sample_seuratobject_20241217.R $path_to_seurat_object $output_dir $so_file_name

echo "Making tmm file"
path_to_seurat_object=$output_dir/$so_file_name
tmm_file_name=igt1_96_gdT_20250109_tmm.Rds
Rscript $SCRIPT_DIR/limma_make_tmm_20241217.R $path_to_seurat_object $output_dir $tmm_file_name

echo "Fitting"
path_to_tmm_object=$output_dir/$tmm_file_name
fit_file_name=igt1_96_gdT_20250109_fit.Rds
Rscript $SCRIPT_DIR/limma_fit_20241217_level2.IGTHT.R $path_to_seurat_object $path_to_tmm_object $output_dir $fit_file_name

echo "OneVsAll"
path_to_fit_object=$output_dir/$fit_file_name
prefix_file_name=OneVsAll
Rscript $SCRIPT_DIR/limma_contrasts_20250121_OneVsAll_notMiniProlif.R $path_to_seurat_object $path_to_tmm_object $path_to_fit_object $output_dir $prefix_file_name

echo "Prolif miniverse"
path_to_fit_object=$output_dir/$fit_file_name
prefix_file_name=Level2Group
Rscript $SCRIPT_DIR/limma_contrasts_20250121_Level2Group.R $path_to_seurat_object $path_to_tmm_object $path_to_fit_object $output_dir $prefix_file_name

