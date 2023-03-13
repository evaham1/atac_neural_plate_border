#!/bin/sh
#SBATCH --job-name=Integrated-npb/NF-cellranger_align
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=eva.hamrud@crick.ac.uk

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.6
ml Singularity/3.6.4
ml Graphviz

export TERM=xterm
export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR=/camp/home/hamrude/briscoej/working/hamrude/NF_singularity

## UPDATE PIPLINE
nextflow pull Streit-lab/cellranger_multiomic

nextflow run Streit-lab/cellranger_multiomic \
    -r test_release \
    -c ./conf/crick_params.config \
    --sample_sheet ./samplesheet.csv \
    --outdir ../output/NF-cellranger_align \
    --email eva.hamrud@crick.ac.uk \
    -resume