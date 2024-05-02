#!/bin/bash

# Specify inputs
gtf_file="$1"
index_file="$2"
output_file="$3"

awk -F'\t' '$3 ~ /gene/' "$gtf_file" > temp.gtf

while IFS=$'\t' read -r -a line; do
    # Extract the chromosome and store it as a variable
    chrom="${line[0]}"

    # Extend 3kb from the gene start site depending on strand
    if [ "${line[6]}" == "+" ]; then
        start=$(( ${line[3]} - 2000 ))
        end="${line[3]}"
    elif [ "${line[6]}" == "-" ]; then
        start="${line[4]}"
        end=$(( ${line[4]} + 2000 ))
    fi

    # Extract the chromosome length from the index file
    chrom_length=$(awk -v chrom="$chrom" '$1 == chrom {print $2}' "$index_file")

    # Use the chromosome length to adjust the promoter coordinates
    if (( start < 0 )); then
        start=0
    fi
    if (( end > chrom_length )); then
        end=$chrom_length
    fi

    # Extract gene_id and gene_name values
    a=()
    IFS='\ ' read -ra a <<< "${line[8]}"

    gene_id=$(echo "${a[1]}" | sed 's/"//g; s/;//g')

    if [[ "${a[4]}" == 'gene_name' ]]; then
        gene_name=$(echo "${a[5]}" | sed 's/"//g; s/;//g')
    else
        gene_name=$gene_id
    fi

    # Print out the elements of interest
    echo -e "chr$chrom\t$start\t$end\t${line[5]}\t${line[6]}\t$gene_id\t$gene_name"
done < temp.gtf > "$output_file"

# Remove the temp.gtf file
rm temp.gtf