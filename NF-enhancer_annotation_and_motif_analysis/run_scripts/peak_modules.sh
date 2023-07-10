#!/bin/sh
#SBATCH --job-name=NF-enhancer_annotation_and_motif_analysis
#SBATCH -t 42:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=hamrude@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Java/11.0.2
ml Nextflow/22.10.3
ml Singularity/3.6.4

export NXF_VER=22.10.3

## UPDATE PIPLINE
nextflow pull Streit-lab/enhancer_annotation_and_motif_analysis

nextflow run Streit-lab/enhancer_annotation_and_motif_analysis \
    -r main \
    -c ./conf/crick_params.config \
    -profile singularity \
    --peaks_bed  ./data/samplesheet.csv \
    --outdir ../output/NF-enhancer_annotation_and_motif_analysis/peak_modules \
    --email hamrude@crick.ac.uk \
    -resume
