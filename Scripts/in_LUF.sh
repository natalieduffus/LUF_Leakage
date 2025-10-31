#!/bin/bash

#SBATCH --job-name=LUF # the name for the cluster scheduler
#SBATCH --time=01:00:00 # Maximum allowed runtime per iteration
#SBATCH --mem-per-cpu=50G
#SBATCH --array=1-1000 # the number of iterations
#SBATCH --output=logfiles/LUF_%A_%a.out # the name of the output files
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.duffus@biology.ox.ac.uk

# load the anaconda module
module load Anaconda3

# Make a string variable corresponding to your $DATA directory, then conda_envs/nats_modelling. E.g. for Dave the absolute path is /data/zool-mosquito_ecology/zool2291/conda_envs/nats_modelling

export CONPREFIX=$DATA/conda_envs/nats_modelling

# activate your new conda environment
source activate $CONPREFIX

Rscript Cluster_Code.R ${SLURM_ARRAY_TASK_ID}