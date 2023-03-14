#!/bin/sh
#SBATCH --job-name=NF-hichip_fetchngs
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=thierya@crick.ac.uk

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/22.10.3
ml Singularity/3.6.4
ml Graphviz

export TERM=xterm
export NXF_VER=22.10.3
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/briscoej/working/thierya/singularity
export NXF_HOME=/nemo/project/home/thierya/working/nextflow

## UPDATE PIPLINE
nextflow pull nf-core/fetchngs

nextflow run nf-core/fetchngs \
    -r 1.9 \
    -c ./conf/crick_params.config \
    --input SRR_Acc_List.txt \
    --outdir ../output/NF-hichip_fetchngs \
    --nf_core_pipeline atacseq \
    --email thierya@crick.ac.uk \
    -resume
