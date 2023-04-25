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
## sbatch split_and_edit_validpairs.sh input_file.allValidPairs output_file.txt 12

# Get the input file name from the command line argument
input_file="$1"

# Get the output file name from the command line argument
output_file="$2"

# Get the number of columns from the command line argument
num_cols="$3"

# Rename input file to .txt
cp "$input_file" input.txt

# Split the input file into smaller chunks
split -l 1000000 input.txt input_part

# Edit the second and fifth column of each chunk to add 'chr' to the chromosome name
for file in input_part*; do
    awk -v num_cols="$num_cols" 'BEGIN{FS=OFS="\t"}{if(NF==num_cols){$2="chr"$2;$5="chr"$5}; print}' "${file}" > "${file}.edited" &
done

# Wait for all editing jobs to finish
wait

# Concatenate the edited chunks into a single output file
cat input_part*.edited > "$output_file"

# Remove the intermediate files
rm input.txt input_part*

# Convert the output file to tab-delimited format
sed -i 's/ /\t/g' "$output_file"

# Print a completion message
echo "Finished editing file."
