#!/bin/sh
#SBATCH --job-name=NF-cellranger_atac_align
#SBATCH -t 2:00:00
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
export NXF_HOME=/nemo/project/home/thierya/.nextflow

## UPDATE PIPLINE
nextflow pull Streit-lab/cellranger_multiomic

nextflow run Streit-lab/cellranger_multiomic \
    -r test_release \
    -c ./conf/crick_params_alex.config \
    --sample_sheet ./samplesheet_alex.csv \
    --outdir ../output/NF-cellranger_align \
    --email thierya@crick.ac.uk \
    -resume