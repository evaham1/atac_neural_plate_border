#!/bin/bash
#SBATCH --job-name=edit_ValidPairs
#SBATCH --output=edit_ValidPairs.out
#SBATCH --error=edit_ValidPairs.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=hamrude@crick.ac.uk

## Bash script to take big .ValidPairs output that comes from HiC-Pro, cut it into smaller pieces, and rename chromosomes from '1' -> 'chr1' etc.
## Usage: sbatch edit_file.sh input_file.allValidPairs output_file.txt 12

# Rename input file to .txt if it has .allValidPairs extension
input_file=$1
if [[ "$input_file" == *.allValidPairs ]]; then
    mv "$input_file" "${input_file/.allValidPairs/.txt}"
    input_file="${input_file/.allValidPairs/.txt}"
fi

# Function to split the input file into smaller chunks
split_file() {
    input_file=$1
    chunk_size=$2
    output_dir=$3
    mkdir -p $output_dir
    split -a 4 -d -l $chunk_size $input_file $output_dir/
}

# Function to edit the second column of a tab-delimited file
edit_file() {
    input_file=$1
    output_file=$2
    num_columns=$3
    awk -F '\t' -v OFS='\t' '{ $2="chr"$2; print }' $input_file > $output_file
}

# Parse command-line arguments
output_file=$2
num_columns=$3

# Split the input file into smaller chunks
chunk_size=1000000 # adjust as needed
output_dir=$(mktemp -d)
split_file $input_file $chunk_size $output_dir

# Edit the second column of each chunk
for chunk_file in $output_dir/*
do
    chunk_output_file=$(basename $chunk_file)
    edit_file $chunk_file $chunk_output_file $num_columns
done

# Combine the edited chunks into a single output file
cat $output_dir/* > $output_file

# Remove the temporary directory
rm -r $output_dir
