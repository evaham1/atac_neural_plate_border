#!/bin/sh
#SBATCH --job-name=NF-hichip_fetchngs
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=thierya@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Java/11.0.2
ml Nextflow/22.10.3
ml Singularity/3.6.4
ml Graphviz

export NXF_VER=22.10.3

## UPDATE PIPLINE
nextflow pull nf-core/fetchngs

nextflow run nf-core/fetchngs \
    -r 1.9 \
    -c ./conf/crick_params.config \
    --input ./data/SRR_Acc_List.txt \
    --outdir ../output/NF-hichip_fetchngs \
    --nf_core_pipeline rnaseq \
    --email hamrude@crick.ac.uk \
    -resume