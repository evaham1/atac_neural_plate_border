#!/bin/bash
#SBATCH --job-name=edit_file
#SBATCH --output=edit_file.out
#SBATCH --error=edit_file.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=hamrude@crick.ac.uk

# call the script to split and edit the file
./split_and_edit_validpairs.sh /flask/scratch/briscoej/hamrude/atac_neural_plate_border/output/NF-hichip/NF-hic/NF_HiChip_r1.txt /flask/scratch/briscoej/hamrude/atac_neural_plate_border/output/NF-hichip/NF-hic/NF_HiChip_r1_edited.txt 12