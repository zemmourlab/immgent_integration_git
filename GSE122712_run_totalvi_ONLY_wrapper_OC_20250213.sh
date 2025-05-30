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

module load gcc/9.2.0 R/4.1.2 python/3.9.14 hdf5/1.14.0 boost/1.75.0 openmpi/4.1.1 fftw/3.3.10 java/jdk-1.8u112 geos/3.10.2
#module load python
#source $(conda info --base)/etc/profile.d/conda.sh

export NUMBA_CACHE_DIR="/tmp"
#source activate /project/zemmour/david/envs/scvi120_20241008
#source ~/scvi-tools_20241105/bin/activate
source /n/groups/cbdm_lab/odc180/Python/envs/scvi_tools_py3.9_250318/bin/activate

SCRIPT_DIR=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/Scripts/
working_dir=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE122712/
path_to_mudata=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/SCVI_Integration/CD8/GSE122712/MERGED_RNA_GSE122712_query_ImmgenT.h5ad
prefix=Trial_8Cores
prefix_SCANVI=Trial_SCANVI_Combined_Fix
batchkey=IGT
categorical_covariate_keys=IGTHT
corrected_counts=False
denoised_data=False
cd $working_dir
python3 $SCRIPT_DIR/run_totalvi_v2_OC.py --working_dir=$working_dir --path_to_mudata=$path_to_mudata --prefix=$prefix --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi.log 2>totalvi.err

module load python/3.8.12

#export NUMBA_CACHE_DIR="/tmp"
source /home/odc180/.bashrc
source ~/pymde/bin/activate
#source activate ~/pymde/bin/activate
#mkdir mde/

mde_ref_file=/n/groups/cbdm_lab/odc180/ImmgenT_workshop/Query_Integration/CD8AB/mde_incremental_CD8AB.csv
totalvi_integrated_file=Trial_8Cores/latent.csv
python $SCRIPT_DIR/run_mde_incremental_OC.py $working_dir $prefix $mde_ref_file $totalvi_integrated_file >mde.log 2>mde.err

deactivate

#module load gcc/9.2.0 python/3.10.11 miniconda3/23.1.0 cuda/12.1 hdf5/1.14.0 boost/1.75.0 openmpi/4.1.1 fftw/3.3.10 java/jdk-1.8u112 geos/3.10.2
#export NUMBA_CACHE_DIR="/tmp"
##source activate /project/zemmour/david/envs/scvi120_20241008
##source ~/scvi-tools_20241105/bin/activate
#source activate scarches

#python $SCRIPT_DIR/run_totalvi_v2_OC_withSCANVI.py --working_dir=$working_dir --path_to_mudata=$path_to_mudata --prefix=$prefix --prefix_SCANVI=$prefix_SCANVI --batchkey=$batchkey --categorical_covariate_keys=$categorical_covariate_keys --corrected_counts=$corrected_counts --denoised_data=$denoised_data >totalvi_withSCANVI.log 2>totalvi_withSCANVI.err

