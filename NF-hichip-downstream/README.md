# NF-hichip-downstream

This folder contains all the code used to perform downstream analysis of publicly avaliable HiChip data from the chick ectoderm ([Azambuja and Simoes-Costa 2021](https://pubmed.ncbi.nlm.nih.gov/33852891/)). The code is organised into a custom Nextflow pipeline. 

For this project, loop calling was performed on the HiChip data to predict enhancer-target-gene pairings (ETG pairs) on 3 samples of HiChip data collected from the chick ectoderm. This was performed using the package [HiCDCPlus](https://github.com/mervesa/HiCDCPlus). 

First, the chick genome was split into equal sized 5kb bins and these bins were overlapped with promoters and predicted enhancers. Second, loop calling was performed between these bins to identify which bins (i.e. regions of the genome) are in close physical proximity. As there were three samples, loops were called for all three indpendently and then all loops were considered together. Third, these loops were explored to identify how many of them represent potential ETG pairs and if these pairs can explain the regulation of key ectodermal genes. 

## 1. Bin generation and overlapping with genes and enhancers
This part of the pipeline is independent of the HiChip data, instead it is dependent on the genome (in this case Galgal6 genome was used) and peak coordinates which are used as a proxy for active enhancers (in this case the coordinates of peaks called from the scATAC-seq chick ectoderm atlas was used). These files can be changed as indicated by the orange colour in the workflow schematic. Once edited, this workflow will emit genomic bin coordinates with information about which bins contain enhancers (i.e. peaks) and genes (i.e. promoters). 

*add schematic here

### Re-using the steps in this workflow
The processes [EXTRACT_PROMOTERS](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/modules/local/extract_promoters/main.nf) and [INTERSECT_BINS](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/modules/local/intersect_bins/main.nf) are custom-made Nextflow modules. The EXTRACT_PROMOTERS module calls a Bash script called [extract_promoters.sh](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/extract_promoters.sh) which identifies each gene start site and extracts 3kb upstream of it as a promoter whilst accounting for different chromsome lengths. The INTERSECT_BINS module runs 'bedtools intersect' from the [bedtools toolset](https://bedtools.readthedocs.io/en/latest/). 

The process GENERATE_BINS is a generic R module, which runs the R script [generate_bins.R](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/generate_bins.R) which generates bins using the Gallus gallus genome. 
 
Each of these processes can be re-used as they are in your own pipeline or edited for different purposes. 

## 2. Loop calling on HiChip data using HiCDCPlus
This part of the pipeline takes the HiChip data which has been aligned to the genome using the NF-core HiC pipeline (see the NF-hic folder in this repository, files end with '.allValidPairs') and calls loops in each HiChip data sample. For this analysis, the loop calling was performed using the chick UCSC genome build, which names each chromosomes as 'chr1', 'chr2', etc. As the ensembl has different genome nomenclature, an extra process was run to edit the HiChip data so it can be processed with the UCSC genome. 

*Picture of workflow

### Re-using the steps in this workflow
The process [EDIT_VALIDPAIR](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/modules/local/edit_ValidPairs/main.nf) is a custom Nextflow module which runs the bash script [edit_validpairs.sh](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/edit_validpairs.sh) which adds 'chr' to the beginning of each chromosome in the data. This module and its corresponding script can be reused into any Nextflow pipeline to add 'chr' to validpairs datasets. 

The process LOOP_CALL is a generic R module, which runs the R script [loop_calling_hicdc.R](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/loop_calling_hicdc.R) which takes the input valid pairs dataset and uses the HiCDCPlus package to call loops. This script can be edited to change the parameters and run on any HiChip dataset. 

## 3. Explore the ETG pairs
The final process, DIFF_LOOPS runs the R script [diff_loops_hicdc.R](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/diff_loops_hicdc.R) to identify the union interactions of all three samples. This process takes a channel that has all three samples combined together. 

### Re-using the steps in this workflow
DIFF_LOOPS runs the R script [diff_loops_hicdc.R](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/diff_loops_hicdc.R) to identify the union interactions of all three samples, however it can be edited as per the HiCDCPlus tutorial to run differential loop calling instead. 

Another R script, [investigate_loops.R](https://github.com/evaham1/atac_neural_plate_border/blob/main/NF-hichip-downstream/bin/investigate_loops.R), has been created to explore the loops and look at their overlap with enhancers the genes. Although this script is not incorporated into the pipeline it can be added using the R module and edited. 