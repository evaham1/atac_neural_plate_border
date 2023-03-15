#!/bin/sh
#SBATCH --job-name=NF-hichip_alignment
#SBATCH -t 42:00:00
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
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/singularity
export NXF_HOME=/nemo/project/home/thierya/.nextflow

## UPDATE PIPLINE
nextflow pull nf-core/hic

nextflow run nf-core/hic \
    -r 2.0.0 \
    -c ./conf/crick_params.config \
    --digestion 'mboi' \
    --input  ./data/samplesheet.csv \
    --outdir ../output/NF-hichip_alignment \
    --email thierya@crick.ac.uk \
    -resume
