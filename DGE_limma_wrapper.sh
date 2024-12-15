#!/bin/bash
#SBATCH -p highmem
#SBATCH -t 4-23:59:59
#SBATCH --mem 600G
#SBATCH -c 4
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH -w, --nodelist=compute-h-17-55
#SBATCH --mail-type=END,FAIL
# Run as sbatch DGE_limma_wrapper.sh /n/groups/cbdm_lab/odc180/ImmgenT_workshop/DGE_limma/CD4 .
# from Broad Run directory

### Load Relevant Modules
module load gcc/9.2.0 R/4.1.2 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7 java/jdk-1.8u112 geos/3.10.2
export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"


### Change to target directory
cd $2

Rscript $1/DGE_limma_Only_CD4.R 2 > limma.err > limma.log

