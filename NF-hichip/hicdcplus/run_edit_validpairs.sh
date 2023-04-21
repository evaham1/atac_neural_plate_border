#!/bin/bash

### because the ValidPairs file is so big need to split it into smaller pieces

# set the number of lines per file
lines_per_file=1000000

# split the input file into smaller files
split -l $lines_per_file "$1" input_file_part_

# count the number of split files
num_files=$(ls -1 input_file_part_* | wc -l)

# run the original script on each file in parallel
for (( i=0; i<$num_files; i++ ))
do
    # construct the input and output file names
    input_file="input_file_part_$(printf "%02d" $i)"
    output_file="output_file_part_$(printf "%02d" $i)"

    # run the original script on the current file in the background
    ./edit_validpairs.sh "$input_file" "$output_file" "$2" &

    # display a progress bar
    echo -ne "Processing file $((i+1)) of $num_files [$input_file] : [\033[1;32m"
    for (( j=0; j<=40; j++ ))
    do
        if [ $((j*100/40)) -le $((i*100/num_files)) ]; then
            echo -ne "\033[1;32m#\033[0m"
        else
            echo -ne "\033[1;31m-\033[0m"
        fi
    done
    echo -ne "] $((i*100/num_files))% \r"
done

# wait for all background processes to finish
wait

# concatenate the output files
cat output_file_part_* > "$3"

# clean up the split files
rm -f input_file_part_* output_file_part_*
