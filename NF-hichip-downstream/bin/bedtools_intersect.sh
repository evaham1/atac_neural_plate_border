#!/bin/bash

# Usage: ./betools_intersect.sh <bin.bed> <feature.bed> <output.bed>

# Get the bin.bed file name from the command line argument
bins="$1"

# Get the feature.bed file name from the command line argument
features="$2"

# Get the output file name from the command line argument
output="$3"

bedtools intersect -a $bins -b $features -wa -wb > $output