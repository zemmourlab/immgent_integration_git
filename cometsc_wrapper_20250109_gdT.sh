#!/bin/bash
#SBATCH -A pi-zemmour ##SBATCH -q jfkfloor2 --exclusive 
#SBATCH --partition=caslake 
#SBATCH --nodes=1 
#SBATCH --mem=96GB
#SBATCH -J cometsc_%j            
#SBATCH -o cometsc_%j.log
#SBATCH -t 24:00:00              ##SBATCH --mem=8G
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL ; 
#SBATCH --mail-user=zemmour@rcc.uchicago.edu   # Email to which notifications will be sent

#run as: sbatch $SCRIPT_DIR/cometsc_wrapper_20250109_gdT.sh

module load python/anaconda-2022.05 
source $(conda info --base)/etc/profile.d/conda.sh
module load hdf5/1.12.0 #for BPCells but also others
module load openblas/0.3.13
source activate /project/zemmour/david/envs/cometsc #torch worjks but scvi not, jax error
module load openblas/0.3.13 #load again or error


SCRIPT_DIR=/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git
path_to_wd=/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/gdT/cometsc/20250109

cd $path_to_wd

Comet input_data_log.txt input_umap.txt input_clusters.txt k2_log -C 36 -K 2

Comet input_data.txt input_umap.txt input_clusters.txt k2 -C 36 -K 2 # -C is for cores, -K is the number of gene combinations for two does single and 2 markers -->

