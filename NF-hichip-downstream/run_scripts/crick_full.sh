#!/bin/bash
#SBATCH --job-name=atac-NPB
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=eva.hamrud@crick.ac.uk

export TERM=xterm
export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR=/nemo/lab/briscoej/working/hamrude/NF_singularity
export NXF_HOME=/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream
export NXF_WORK=work/

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.6
ml Singularity/3.4.2
ml Graphviz

##nextflow run ./main.nf -dump-hashes \
nextflow run ./main.nf \
--outdir ../output/NF-hichip-downstream \
-profile crick_full \
-resume