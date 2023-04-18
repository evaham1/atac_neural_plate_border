#!/bin/sh
#SBATCH --job-name=NF-FitHiChiP
#SBATCH -t 42:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=hamrude@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Java/11.0.2
ml Singularity/3.6.4

## RUN FITHICHIP BASH SCRIPT
bash ../run_fithichip/FitHiChIP/FitHiChIP_Singularity.sh -C ../run_fithichip/data/configfile_NF_r1
bash ../run_fithichip/FitHiChIP/FitHiChIP_Singularity.sh -C ../run_fithichip/data/configfile_NF_r2
