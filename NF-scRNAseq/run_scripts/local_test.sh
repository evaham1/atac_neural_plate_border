#!/bin/bash

nextflow run ./main.nf \
--outdir ../output/NF-scRNAseq \
-profile local_test \
-resume