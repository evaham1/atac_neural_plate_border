#!/bin/bash

# Specify inputs
gtf_file="$1"
output_file="$2"

# Extracts all genes
# Identifies orientation of gene and then extract promoter which is gene start +/- 2kb
# Changes chromosome names to 'chr1'
# Extracts gene name or if there is no gene name, gene ID
awk -F'\t' '$3 ~ /gene/ {
    if ($7 == "+") {
        $5 = $4
        $4 = $4 - 2000
    } else if ($7 == "-") {
        $4 = $5
        $5 = $5 + 2000
    }
    # Extract gene_id and gene_name values within quotations
    split($9, a, "\"");
    gene_id=a[2];
         if (a[10]=="") { 
             gene_name=gene_id;
         } else { 
             gene_name=a[6];
         }
    # Create output with selected columns and added gene_id and gene_name variables
    output = "chr"$1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" gene_id "\t" gene_name
    print output
}' "$gtf_file" > "$output_file"
rm "$gtf_file"