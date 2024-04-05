#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@fht.org
#SBATCH --job-name pqtl_pipeline
#SBATCH --output %j_pqtl_pipeline.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 8G
#SBATCH --time 48:00:00

source ~/.bashrc
module -s load singularity/3.8.5

# set some singularity directories depending on frontend/computing node/vm
case $(hostname) in
  hnode*)
    export SINGULARITY_TMPDIR=/tmp/
    export SINGULARITY_BIND="/cm,/exchange,/processing_data,/project,/scratch,/center,/group,/facility,/ssu"
    ;;
  cnode*|gnode*)
    export SINGULARITY_TMPDIR=$TMPDIR
    export SINGULARITY_BIND="/cm,/exchange,/processing_data,/project,/localscratch,/scratch,/center,/group,/facility,/ssu"
    ;;
  lin-hds-*)
    export SINGULARITY_TMPDIR=/tmp/
    export SINGULARITY_BIND="/processing_data,/project,/center,/group,/facility,/ssu,/exchange"
    ;;
  *)
    export SINGULARITY_TMPDIR=/var/tmp/
    ;;
esac

# change to a workspace in scratch
cd /scratch/$USER
mkdir -p pqtl_pipeline
cd pqtl_pipeline

# set a project name and clone the pipeline
PROJECT_NAME="test"
if [ ! -d ${PROJECT_NAME} ]; then
  git clone --recurse-submodules https://github.com/ht-diva/pqtl_pipeline.git ${PROJECT_NAME}
fi
cd ${PROJECT_NAME}

# copy reference data
mkdir -p data
rsync -avz /exchange/healthds/pQTL/pQTL_workplace/TEST_INTERVAL/BATCH/INTERVAL_NonImp_residuals_final.txt data/
rsync -avz /exchange/healthds/public_data/gwaslab_reference_dataset/* data/
rsync -avz /exchange/healthds/public_data/reference_genomes/GRCh37/GCA_000001405.14_GRCh37.p13_full_analysis_set.fna data/

# run the pipeline
conda activate /exchange/healthds/software/envs/snakemake
make run
