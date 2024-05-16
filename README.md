# Neural Plate Border Analysis

This repository contains all the code used for analysis pertaining to the Neural Plate Border project. This project focused on understanding how cell fate decisions are made in the developing ectoderm of the chick embryo. This project involved analysis of scRNA-seq and scATAC-seq data created by the Streit lab, and publicly avaliable CUT&RUN and HiChip data. 

## Repository Organisation
Each folder beginning with 'NF-' contains an independent subset of the analysis. Some of the folders run NF-core pipelines (click [here](https://nf-co.re/pipelines) to see all NF-core pipelines), and some of the folders run our shared Streitlab Nextflow pipelines (click [here](https://github.com/Streit-lab) to see all Streitlab repositories). Two folders (NF-hichip_downstream and NF-downstream_analysis) contain custom Nextflow pipelines. 

| Folder | Pipeline Origin | Description | README | 
| ------ | ------ | ------ | ------ | 
| NF-cellranger_align | Streit-lab pipeline | Aligns scRNA-seq and scATAC-seq data  | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cellranger_align) | 
| NF-cut_and_run | NF-core pipeline | Downloads and aligns CUT&RUN data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-cutandrun) | 
| NF-downstream_analysis | Custom pipeline | Analyses scATAC-seq data and integrate it with previously published scRNA-seq data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-downstream_analysis) | 
| NF-enhancer_annotation_and_motif_analysis | Streit-lab pipeline | Annotates genomic regions to genes and run motif scanning | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-enhancer_annotation_and_motif_analysis) | 
| NF-hic | NF-core pipeline | Aligns HiChip daata | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hic)
| NF-hichip_downstream | Custom pipeline | Analyses HiChip data | [README file](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream) | 

## Reproducing the analysis of this project
As stated above, different subsets of this project's analysis have been run in different ways. 

To re-run the analysis that utilise the Streit-lab pipelines or NF-core pipelines please see the README files linked in the table above and consult the documentation of the original pipelines. 

To reproduce the analysis performed in the two custom pipelines: NF-downstream_analysis and NF-hichip_downstream, perform the following steps:
1) Clone this repository from Github
2) Ensure Nextflow is installed on the computer or HPC in which the pipeline will be run (HPC or cloud server recommended as computing requirements are high)
3) Download and process the raw data (using the Nextflow hic pipeline for the HiChip or Cellranger for the scATAC-seq data) and change the samplesheets to reflect the correct path to the aligned data. 
[Samplesheet for NF-hichip_downstream input](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/samplesheet_validpairs.csv) / [Samplesheet for NF-downstream_analysis input](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/samplesheets/samplesheet_aligned.csv)
4) The bash scripts in the 'run_scripts' folder are how the Nextflow pipelines are triggered. These may need to be edited depending on how modules are loaded in the HPC and to set cache directories for singularity containers and nextflow. 
[Run script for NF-hichip_downstream pipeline](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/run_scripts/crick_full.sh) / [Run script for NF-downstream_analysis pipeline](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/run_scripts/crick_full.sh)
5) The run script also specifies a 'profile' which includes platform-specific parameters and paths to reference genomes. A new profile can be created, named and specified in the run script. This [example profile](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/conf/crick_test.config) can be used as a template.
6) Once the pipeline has been adapted to run from a different computer or HPC system, it can be executed by running the run script

## Reusing code from this repository for other projects
The two custom Nextflow pipelines, NF-downstream_analysis and NF-hichip_downstream, can be reused in whole or in part for other projects and analyses. Below is instructions on how the whole pipeline, modules of the pipeline or scripts of the pipeline can be re-used and re-purposed. For more information on what the different steps of these pipelines do and which might be useful for re-purposing see their README docs ([NF-downstream_analysis README](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-downstream_analysis), [NF-hichip_downstream README](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream)). 

### Reusing whole Nextflow pipeline with different input data - NF_hichip-downstream
Despite being made of some generalised components, the whole NF-downstream_analysis pipeline is specific to the data and question of the NPB project, and so it is not recommeded to run these pipelines entirely on different data. 

In contrast, the NF-hichip_downstream pipeline can be used to predict enhancer-gene pairs in any HiChip dataset. To run the pipeline with different data, 3 data inputs need to be changed. 
1) the [samplesheet](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/samplesheet_validpairs.csv) should be modified to reflect the different HiChip data in the form of Valid Pairs - this data must have been pre-processed using the NF-core HiC pipeline. 
2) The genome should be adapted to the species of interest by changing the path to the gtf file and genome index file in the [profile config](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/conf/crick_full.config). 
3) The 'peaks' path should be a bed file of where ehancers are predicted to be located. In this project the peaks from the scATAC-seq data were used, however bulk peaks, Chip-seq peaks or equivalent can be used. This path is also specified in the [profile config](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/conf/crick_full.config).

For more information on what the NF-hichip_downstream pipeline does see the [pipeline's documentation](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream).

### Reusing Nextflow modules
Nextflow [modules](https://www.nextflow.io/docs/latest/module.html) are called into the pipeline as processes. The modules used for the NF_hichip-downstream pipeline can be found [here](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream/modules/local) and for the NF-downstream_analysis pipeline [here](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-downstream_analysis/modules/local). These modules can be re-used in other Nextflow pipelines by copying their main.nf file and loading them into a pipeline as a process 
e.g.
```
include {EXTRACT_PROMOTERS} from "$baseDir/modules/local/extract_promoters/main"
``` 
Looking at the main.nf Nextflow script for the module can reveal what the module does and what it expects as input and output, for example see the [module to extract promoters from a gtf](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/modules/local/extract_promoters/main.nf). 

Two of the Nextflow modules used in these pipelines are generic in that they run either an [R script](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/modules/local/r/main.nf) or a [python script](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-downstream_analysis/modules/local/python/main.nf). Please note that both modules expect a tuple as input data - i.e. the data and associated metadata. 

> Note that each Nextflow module is associated with a Docker container,
> this can be seen in the module's main.nf script. These containers 
> overcome dependency issues and make the Nextflow pipeline portable. 
> Some of the containers are from the Docker repository Docker.hub
> whilst others are custom-made, and their Docker scripts can be 
> found in the same folder as the module's main.nf script

### Reusing R and python scripts
The generic R and python modules mentioned above can be read into the Nextflow pipeline together with a specific script like so:
```
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_peak_calling.R", checkIfExists: true) )
```

The custom pipelines in this repository include many R and python scripts that can be found in the 'bin' folder of each pipeline. Most of these scripts are highly generalisable and can be used to analyse any scATAC-seq, scRNA-seq and HiChip data. Most of the scripts contain a list of script options at the top which can be overwritten to adjust parameters or skip parts of the script, these options look like this:
```
# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("", "--stage_clust_res"), action = "store", type = "double", help = "clustering resolution for stage data", default = 1),
    make_option(c("", "--full_clust_res"), action = "store", type = "double", help = "clustering resolution for full data", default = 2),
    make_option(c("", "--clustree_stage"), action = "store", type = "logical", help = "whether to run clustree plot on stage data", default = FALSE),
    make_option(c("", "--clustree_full"), action = "store", type = "logical", help = "whether to run clustree plot on full data", default = FALSE),
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)
```
and can be overwritten in a config file if running the script using the R/python modules in Nextflow pipeline, or they can just be added as arguments if the script is being run from the commandline. 

For more information on each R and python script and what they do see the pipeline's README docs ([NF-downstream_analysis README](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-downstream_analysis), [NF-hichip_downstream README](https://github.com/evaham1/atac_neural_plate_border/tree/main/NF-hichip-downstream)). 
