#!/bin/bash
#SBATCH -p medium
#SBATCH -t 2-12:00:00
#SBATCH --mem 64G
#SBATCH -c 8
#SBATCH -o wrapper_totalvi.log
#SBATCH -e wrapper_totalvi.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch run_totalvi_wrapper_OC_20250213.sh .
# from working directory

module load gcc/9.2.0 python/3.9.14 R/4.1.2
#module load python
#source $(conda info --base)/etc/profile.d/conda.sh

python3 -c "print('Hello from Python')"

#export NUMBA_CACHE_DIR="/tmp"
#source activate /project/zemmour/david/envs/scvi120_20241008
#source ~/scvi-tools_20241105/bin/activate
source /n/groups/cbdm_lab/odc180/Python/envs/scvi_tools_pymde_annoy_py3.9_20250522/bin/activate

SCRIPT_DIR=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/
working_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE122712/Webpage_Trial/
path_to_ImmgenT=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Webpage_Trial/ImmgenT_downsampled_20250505_adata.h5ad
path_to_query=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Webpage_Trial/query_adata.h5ad
path_to_anndata=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE122712/MERGED_RNA_GSE122712_query_ImmgenT.h5ad
prefix=Full_Trial
batchkey=IGT
categorical_covariate_keys=IGTHT
corrected_counts=False
denoised_data=False
mde_ref_file=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/ImmgenT_freeze_20250109/igt1_104_withtotalvi20250505_downsampled_level1_and_level2_mde.csv
cd $working_dir
/n/groups/cbdm_lab/odc180/Python/envs/scvi_tools_pymde_annoy_py3.9_20250522/bin/python3 -c "print('Hello from System Python')"

python3 $SCRIPT_DIR/webpage_integration_anndata_SCANVI_pymde_script.py --working_dir=$working_dir --path_to_ImmgenT=$path_to_ImmgenT --path_to_query=$path_to_query --path_to_anndata=$path_to_anndata --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data --mde_ref_file=$mde_ref_file >totalvi.log 2>totalvi.err

export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
output_file=$working_dir/$prefix/output_file.csv
user_output_file=$working_dir/$prefix/user_output_file.csv
Rscript $SCRIPT_DIR/Plotting_for_webpage.R $output_file $user_output_file $prefix >plot.log 2>plot.err

