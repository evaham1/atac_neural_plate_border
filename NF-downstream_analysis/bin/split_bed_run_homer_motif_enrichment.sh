#!/bin/bash

# Input file provided as an argument
input_file="$1"

# Initialize an array to store the unique PM prefixes
unique_prefixes=()

# Read the input file line by line
while IFS=$'\t' read -r col1 col2 col3 pm_value; do
  # Extract the prefix before the hyphen in the PM value
  pm_prefix="${pm_value%%-*}"
  
  # Determine the output file based on the PM prefix
  output_file="${pm_prefix}.bed"
  
  # Check if the output file already exists in the array
  if ! [[ " ${unique_prefixes[*]} " =~ " $pm_prefix " ]]; then
    # Initialize the output file if it doesn't exist
    > "$output_file"
    # Add the prefix to the array to mark it as handled
    unique_prefixes+=("$pm_prefix")
  fi
  
  # Write the data line to the corresponding output file
  echo -e "$col1\t$col2\t$col3\t$pm_value" >> "$output_file"

  # Run homer with each of the bed files
  findMotifsGenome.pl "$output_file" "$2" "$pm_prefix" -size given

done < "$input_file"