#!/bin/sh
#SBATCH --job-name=NF-enhancer_annotation_and_motif_analysis
#SBATCH -t 42:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=hamrude@crick.ac.uk

export TERM=xterm
export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/briscoej/home/users/hamrude/singularity
export NXF_HOME=/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-enhancer_annotation_and_motif_analysis
export NXF_WORK=work/

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
    -profile test,singularity \
    --fasta /nemo/lab/briscoej/home/users/hamrude/raw_data/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
    --gtf /nemo/lab/briscoej/home/users/hamrude/raw_data/genomes/galgal6/tag_chroms.gtf \
    --outdir output