#!/bin/bash

# Specify inputs
input_file="$1"
output_file="$2"

# Rename input file to .txt
mv "$input_file" "${input_file%.*}.txt"

# Split the input file into smaller chunks
split -l 1000000 "${input_file%.*}.txt" input_part

# Edit the second and fifth column of each chunk to add 'chr' to the chromosome name
for file in input_part*; do
    awk 'BEGIN{FS=OFS="\t"}{$2="chr"$2;$5="chr"$5}1' "${file}" > "${file}.edited" &
done

# Wait for all editing jobs to finish
wait

# Concatenate the edited chunks into a single output file
cat input_part*.edited > "$output_file"

# Remove the intermediate files
rm input_part* "${input_file%.*}.txt"

# Remove the input folder
rm -r input

# Convert the output file to tab-delimited format and rename
sed -i 's/ /\t/g' "$output_file"