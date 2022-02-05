#!/bin/bash
#SBATCH --job-name=atac-NPB
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=eva.hamrud@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.3
ml Singularity/3.4.2
ml Graphviz

export NXF_VER=21.10.3

nextflow run ./main.nf \
--outdir ../output/NF-downstream_analysis \
-profile crick_full \
-resume