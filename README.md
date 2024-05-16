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

### Reproducing entire pipeline without modification
To reproduce the entire analysis of this project without modification, only two things need to be changed: the paths to the data need to be modified so they are relative to where you are running the pipeline from, and the profile config needs to be adapted to reflect the computing resources avaliable. 

1) Clone this repository from Github
2) Ensure Nextflow is installed on the computer or HPC in which the pipeline will be run (HPC or cloud server recommended as computing requirements are high)
3) Download and process the raw data (using the Nextflow hic pipeline for the HiChip or Cellranger for the scATAC-seq data) and change the samplesheets to reflect the correct path to the aligned data. 
[Samplesheet for NF-hichip_downstream input](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/samplesheet_validpairs.csv) / [Samplesheet for NF-downstream_analysis input](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/samplesheets/samplesheet_aligned.csv)
4) The bash scripts in the 'run_scripts' folder are how the Nextflow pipelines are triggered. These may need to be edited depending on how modules are loaded in the HPC and to set cache directories for singularity containers and nextflow. 
[Run script for NF-hichip_downstream pipeline](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/run_scripts/crick_full.sh) / [Run script for NF-downstream_analysis pipeline](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/run_scripts/crick_full.sh)
5) The run script also specifies a 'profile' which includes platform-specific parameters and paths to reference genomes. A new profile can be created, named and specified in the run script. This [example profile](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/conf/crick_test.config) can be used as a template.
6) Once the pipeline has been adapted to run from a different computer or HPC system, it can be executed by running the run script

### Reusing whole Nextflow pipeline with different input data - NF_hichip-downstream
Despite being made of some generalised components, the whole NF-downstream_analysis pipeline is specific to the data and question of the NPB project, and so it is not recommeded to run these pipelines entirely on different data. 

The NF-hichip_downstream pipeline can be used to predict enhancer-gene pairs in any HiChip dataset. To run the pipeline with different data, 3 data inputs need to be changed. 
1) the [samplesheet](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/samplesheet_validpairs.csv) should be modified to reflect the different HiChip data in the form of Valid Pairs - this data must have been pre-processed using the NF-core HiC pipeline. 
2) The genome should be adapted to the species of interest by changing the path to the gtf file and genome index file in the [profile config](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/conf/crick_full.config). 
3) The 'peaks' path should be a bed file of where ehancers are predicted to be located. In this project the peaks from the scATAC-seq data were used, however bulk peaks, Chip-seq peaks or equivalent can be used. This path is also specified in the [profile config](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/conf/crick_full.config).

For more information on what the NF-hichip_downstream pipeline does see the [pipeline's documentation](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream).

### Reusing Nextflow modules
Nextflow [modules](https://www.nextflow.io/docs/latest/module.html) are called into the pipeline as processes. The modules used for the NF_hichip-downstream pipeline can be found here and for the NF-downstream_analysis pipeline here. 

### Reusing R and python scripts


how to just reproduce everything exactly - just change run script and profile
how to adjust things like data and params, see read mes for each one
can also just use the r or python scripts independently, have been made as generalisable as possible (eg converting seurat/scanpy scripts)

data is read in using sample sheets for multiple samples and the profile with channels for params/reference genomes/etc
explain how to edit these if just changing input data

