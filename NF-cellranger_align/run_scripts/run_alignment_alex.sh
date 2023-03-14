#!/bin/sh
#SBATCH --job-name=NF-cellranger_atac_align
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=thierya@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/22.10.3
ml Singularity/3.4.2
ml Graphviz

export NXF_VER=22.10.3
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/singularity

## UPDATE PIPLINE
nextflow pull Streit-lab/cellranger_multiomic

nextflow run Streit-lab/cellranger_multiomic \
    -r test_release \
    -c ./conf/crick_params.config \
    --sample_sheet ./samplesheet_alex.csv \
    --outdir ../output/NF-cellranger_align \
    --email eva.hamrud@crick.ac.uk \
    -resume