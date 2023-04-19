#!/bin/sh
#SBATCH --job-name=NF-hichip_alignment
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
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/briscoej/working/hamrude/NF_singularity
export NXF_HOME=/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-cutandrun
export NXF_WORK=work/

## UPDATE PIPLINE
nextflow pull nf-core/cutandrun

nextflow run nf-core/cutandrun \
    -r dev \
    -c ./conf/crick_params.config \
    --input  ./data/samplesheet.csv \
    --normalisation_mode 'CPM' \
    --peakcaller 'macs2' \
    --outdir ../output/NF-cutandrun-H3k27ac \
    --email hamrude@crick.ac.uk \
    -resume
