#!/bin/bash

# Specify inputs
gtf_file="$1"
output_file="$2"

awk -F'\t' '$3 ~ /exon/' "$gtf_file" > temp.gtf

while IFS=$'\t' read -r -a line; do
    # Extract the chromosome, start, end and strand and store it as a variable
    chrom="${line[0]}"
    start="${line[3]}"
    end="${line[4]}"
    strand="${line[6]}"

    # Extract gene_id, transcript_id and gene_name
    a=()
    IFS='\ ' read -ra a <<< "${line[8]}"

    gene_id=$(echo "${a[1]}" | sed 's/"//g; s/;//g')
    transcript_id=$(echo "${a[5]}" | sed 's/"//g; s/;//g')

    if [[ "${a[10]}" == 'gene_name' ]]; then
        gene_name=$(echo "${a[11]}" | sed 's/"//g; s/;//g')
    else
        gene_name=$gene_id
    fi

    # Extract transcript_id

    # Print out the elements of interest
    #echo -e "chr$chrom\t$start\t$end\t$strand\t$transcript_id\t$gene_id\t$gene_name"
    echo -e ""chr"$chrom\t$start\t$end\t$strand\t$gene_id\t$transcript_id\t$gene_name"
done < temp.gtf > "$output_file"

# Remove the temp.gtf file
rm temp.gtf