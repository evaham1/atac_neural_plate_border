# Neural Plate Border Analysis

This repository contains all the code used for analysis pertaining to the Neural Plate Border project. This project focused on understanding how cell fate decisions are made in the developing ectoderm of the chick embryo. This project involved analysis of scRNA-seq and scATAC-seq data created by the Streit lab, and publicly avaliable CUT&RUN and HiChip data. 

## Repository Organisation
Each folder beginning with 'NF-' contains an independent subset of the analysis. Some of the folders run NF-core pipelines (click [here](https://nf-co.re/pipelines) to see all NF-core pipelines), and some of the folders run our shared Streitlab Nextflow pipelines (click [here](https://github.com/Streit-lab) to see all Streitlab repositories). Two folders (NF-hichip_downstream and NF-downstream_analysis) contain custom Nextflow pipelines, details of how to re-run these pipelines can be found below. 

| Folder | Pipeline Origin | Description | README | 
| ------ | ------ | ------ | ------ | 
| NF-cellranger_align | Streit-lab pipeline | Aligns scRNA-seq and scATAC-seq data  | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cellranger_align) | 
| NF-cut_and_run | NF-core pipeline | Downloads and aligns CUT&RUN data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cutandrun) | 
| NF-downstream_analysis | Custom pipeline | Analyses scATAC-seq data and integrate it with previously published scRNA-seq data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-downstream_analysis) | 
| NF-enhancer_annotation_and_motif_analysis | Streit-lab pipeline | Annotates genomic regions to genes and run motif scanning | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-enhancer_annotation_and_motif_analysis) | 
| NF-hic | NF-core pipeline | Aligns HiChip daata | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hic)
| NF-hichip_downstream | Custom pipeline | Analyses HiChip data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream) | 

## Reproducing analysis in NF-hichip_downstream and NF-downstream_analysis
Although the two custom pipelines in this repository have been created with a specific analysis in mind, their organisation in a Nextflow pipeline means that they can be easily reused in whole or in part. 

how to just reproduce everything exactly - just change run script and profile
how to adjust things like data and params, see read mes for each one
can also just use the r or python scripts independently, have been made as generalisable as possible (eg converting seurat/scanpy scripts)