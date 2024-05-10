# Neural Plate Border Analysis

This repository contains all the code used for the analysis of scRNA-seq and scATAC-seq data from the chick neural plate border, as well as publically avaliable CUT&RUN and HiChip data. Each folder beginning with 'NF-' contains an independent Nextflow pipeline for a subset of the analysis. Some of the analysis includes custom pipelines whilst others utilise NF-core or shared pipelines. All of the custom code used to analyse the scATAC-seq neural plate border atlas is found in 'NF-downstream_analysis'. 

| Folder | Description | README | 
| ------ | ------ | ------ | 
| NF-cellranger_align | Runs a shared streit-lab alignment pipeline  | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cellranger_align) | 
| NF-cut_and_run | Runs NF-core pipeline to align CUT&RUN data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cellranger_align) | 
| Google Drive | [plugins/googledrive/README.md][PlGd] | README | 
| OneDrive | [plugins/onedrive/README.md][PlOd] | README | 
| Medium | [plugins/medium/README.md][PlMe] | README | 
| Google Analytics | [plugins/googleanalytics/README.md][PlGa] | README | 