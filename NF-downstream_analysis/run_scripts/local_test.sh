#!/bin/bash

nextflow run ./main.nf \
--outdir ../output/NF-downstream_analysis \
-profile local_test \
-resume