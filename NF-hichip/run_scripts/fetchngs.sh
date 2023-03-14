#!/bin/sh
#SBATCH --job-name=NF-hichip_fetchngs
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=thierya@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.6
ml Singularity/3.4.2
ml Graphviz

export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/singularity

## UPDATE PIPLINE
nextflow pull nf-core/fetchngs

nextflow run nf-core/fetchngs \
    -r 1.8 \
    -c ./conf/crick_params.config \
    --input SRR_Acc_List.txt \
    --outdir ../output/NF-hichip_fetchngs \
    --nf_core_pipeline atacseq \
    --email thierya@crick.ac.uk \
    -resume
