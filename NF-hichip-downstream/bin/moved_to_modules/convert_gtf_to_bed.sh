#!/bin/bash

# Usage: ./gtf_to_bed.sh <gtf_file>

if [ $# -eq 0 ]; then
  echo "Usage: ./gtf_to_bed.sh <gtf_file>"
  exit 1
fi

gtf_file="$1"

# Removes header lines starting with "#"
# Only keeps lines with "gene" in the 3rd column
# Splits the 9th column by ";" and then by "\""
# Prints the 1st, 4th, 5th, 9th, 10th, and 7th columns (Chr, start, end, gene_id, gene_name, strand)
# The gene name is the 6th column if it exists, otherwise it is the gene id
# Chromosome names are prepended with "chr"
sed '/^#/d' "$gtf_file" \
| awk -F "\t" '{
    if ($3=="gene") {
        split($9, a, "\"");
        gene_id=a[2];
        chr="chr"$1
        if (a[10]=="") { 
            gene_name=gene_id;
        } else { 
            gene_name=a[6];
        } 
    print chr"\t"$4-1"\t"$5"\t"gene_id"\t"gene_name"\t"$7;
  }
}' \
> "${gtf_file%.*}.bed"