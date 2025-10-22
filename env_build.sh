#!/bin/bash

#SBATCH --job-name=env_build # the name for the cluster scheduler
#SBATCH --time=00:30:00 # Maximum allowed runtime per iteration
#SBATCH --mem-per-cpu=7G
#SBATCH --output=logfiles/env_build_%A.out # the name of the output files
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.duffus@biology.ox.ac.uk

# load the anaconda module
module load Anaconda3

# Make a string variable corresponding to your $DATA directory, then conda_envs/nats_modelling. E.g. for Dave the absolute path is /data/zool-mosquito_ecology/zool2291/conda_envs/nats_modelling

export CONPREFIX=$DATA/conda_envs/nats_modelling
# create an empty conda environment at this path.
# N.B. for this to work you'll need to have a directory in your $DATA directory
# called conda_envs. Create it with mkdir $DATA/conda_envs if you don't have one yet.

conda create --prefix $CONPREFIX

# activate your new conda environment
source activate $CONPREFIX


# accept terms for the conda channels we need to use to install R and the R packages
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main

# install all of your desired packages, from the channels conda-forge, defaults
# and  r. The --yes flag confirms in advance that you want to install them, as 
# otherwise conda assumes that you're working interactively and asks for you to confirm.

conda install R r-sf r-tibble r-dplyr r-units r-readr -c conda-forge -c defaults -c r --yes
