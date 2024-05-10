# NF-cellranger_align

This folder contains the run script to run the alignment of 10X single cell RNA and ATAC sequencing using a [custom Streit-lab Nextflow pipeline](https://github.com/Streit-lab/cellranger_multiomic). This pipeline runs [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) to align sequenced reads to the genome and identify gene counts (scRNA-seq) or fragments and peaks (for scATAC-seq). 

For your own 10X scRNA-seq alignment, we would recommend the [scRNA-seq Nextflow pipeline](https://nf-co.re/scrnaseq/2.6.0) as it is regularly maintained and updated. 