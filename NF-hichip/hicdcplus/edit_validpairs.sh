#!/bin/bash

# Check if the input file, output file, and number of columns arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_file> <output_file> <num_columns>"
    exit 1
fi

# Set the input and output file names as the first and second arguments, respectively
input_file=$1
output_file=$2

# Set the number of columns in the input file as the third argument
num_columns=$3

# Get the total number of lines in the input file
num_lines=$(wc -l < "${input_file}")

# Initialize the line counter and progress bar
line_count=0
progress_bar=""

# Loop through each line of the input file and append "chr" to the beginning of the second column
while read line; do
    # Split the line into columns using the tab delimiter
    IFS=$'\t' read -ra cols <<< "$line"

    # Check if the number of columns in the line matches the expected number of columns
    if [ ${#cols[@]} -ne $num_columns ]; then
        echo "Error: Invalid number of columns in line '$line'"
        exit 1
    fi

    # Modify the second column by appending "chr" to its beginning
    cols[1]="chr${cols[1]}"

    # Join the columns back into a single tab-delimited line and write it to the output file
    echo -e "${cols[@]}" >> "${output_file}"

    # Update the line counter and progress bar
    line_count=$((line_count + 1))
    progress_bar=$(echo "scale=2; $line_count * 100 / $num_lines" | bc -l)

    # Print the progress bar
    printf "Progress: [%-50s] %d%%\r" "${progress_bar// /=}" "${progress_bar%.*}"
done < "${input_file}"

# Print a newline character after the progress bar
printf "\n"