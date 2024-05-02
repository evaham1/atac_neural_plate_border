#!/bin/bash

# Specify inputs
input_file="$1"
output_file="$2"

# Filter chromosomes to only keep those in chroms 1-33, W, Z and MT
# removes "chrAADN05000525.1" "chrAADN05000623.1" "chrAADN05001159.1"
# [36] "chrAADN05001201.1" "chrAADN05001210.1" "chrAADN05001214.1" "chrAADN05001268.1" "chrAADN05001289.1"
# [41] "chrAADN05001290.1" "chrAADN05001317.1" "chrAADN05001329.1" "chrAADN05001331.1" "chrAADN05001366.1"
# [46] "chrAADN05001383.1" "chrAADN05001388.1" "chrAADN05001426.1" "chrAADN05001455.1" "chrAADN05001464.1"
# [51] "chrAADN05001473.1" "chrKZ626819.1"     "chrKZ626825.1"     "chrKZ626826.1"     "chrKZ626828.1"    
# [56] "chrKZ626830.1"     "chrKZ626831.1"     "chrKZ626832.1"     "chrKZ626833.1"     "chrKZ626834.1"    
# [61] "chrKZ626835.1"     "chrKZ626836.1"     "chrKZ626837.1"     "chrKZ626838.1"     "chrKZ626839.1"
awk 'length($2) <= 2' "$input_file" > temp.txt

# Change chromosome names to chr1, chr2, etc.
awk 'BEGIN{FS=OFS="\t"}{$2="chr"$2;$5="chr"$5}1' temp.txt > "$output_file"

# Convert the output file to tab-delimited format and rename
sed -i 's/ /\t/g' "$output_file"

# Remove the temp file
rm temp.txt

















# # Rename input file to .txt
# mv "$input_file" "${input_file%.*}.txt"

# # Split the input file into smaller chunks
# split -l 1000000 "${input_file%.*}.txt" input_part

# # Edit the second and fifth column of each chunk to add 'chr' to the chromosome name
# for file in input_part*; do
#     awk 'BEGIN{FS=OFS="\t"}{$2="chr"$2;$5="chr"$5}1' "${file}" > "${file}.edited" &
# done

# # Wait for all editing jobs to finish
# wait

# # Concatenate the edited chunks into a single output file
# cat input_part*.edited > "$output_file"

# # Remove the intermediate files
# rm input_part* "${input_file%.*}.txt"

# # Remove the input folder
# rm -r input

# # Convert the output file to tab-delimited format and rename
# sed -i 's/ /\t/g' "$output_file"