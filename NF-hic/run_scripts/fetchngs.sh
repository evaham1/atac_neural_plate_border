#!/bin/sh
#SBATCH --job-name=NF-hichip_fetchngs
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=hamrude@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Java/11.0.2
ml Nextflow/22.10.3
ml Singularity/3.6.4
ml Graphviz

export NXF_VER=22.10.3
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/briscoej/working/hamrude/NF_singularity
export NXF_HOME=/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hic
export NXF_WORK=work/

## UPDATE PIPLINE
nextflow pull nf-core/fetchngs

nextflow run nf-core/fetchngs \
    -r 1.9 \
    -c ./conf/crick_params.config \
    --input ./data/SRR_Acc_List.txt \
    --outdir ../output/NF-hichip_fetchngs \
    --email hamrude@crick.ac.uk \
    --nf_core_pipeline hic \
    -resume
