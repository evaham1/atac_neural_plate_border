# Neural Plate Border Analysis

This repository contains all the code used for analysis pertaining to the Neural Plate Border project. 

## Repository Organisation
This repository contains all the code used for the analysis of scRNA-seq and scATAC-seq data from the chick neural plate border, as well as publically avaliable CUT&RUN and HiChip data. Each folder beginning with 'NF-' contains an independent Nextflow pipeline for a subset of the analysis. Some of the analysis includes custom pipelines whilst others utilise NF-core or shared pipelines. All of the custom code used to analyse the scATAC-seq neural plate border atlas is found in 'NF-downstream_analysis'. 

| Folder | Description | README | 
| ------ | ------ | ------ | 
| NF-cellranger_align | Runs a shared streit-lab alignment pipeline  | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cellranger_align) | 
| NF-cut_and_run | Runs NF-core pipeline to align CUT&RUN data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cutandrun) | 
| NF-downstream_analysis* | Runs a custom Nextflow pipeline to perform analysis on scATAC-seq data and integrate it with previously published scRNA-seq data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-downstream_analysis) | 
| NF-enhancer_annotation_and_motif_analysis | Runs a shared streit-lab pipeline to annotate genomic regions to genes and run motif scanning | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-enhancer_annotation_and_motif_analysis) | 
| NF-hic | Runs NF-core pipeline to align HiC daata | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hic)
| NF-hichip_downstream* | Runs custom Nextflow pipeline to perform analysis on HiChip data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream) | 

Pipelines with * indiciate those which are custom Nextflow pipelines created specifically for analysis pertaining to the Neural Plate Border project. Although these pipelines have been created with a specific analysis in mind, their organisation in a Nextflow pipeline means that they can be easily reused in whole or in part (see below). 

## Reproducing analysis
