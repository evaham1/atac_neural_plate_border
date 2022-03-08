#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {EDIT_GTF} from "$baseDir/modules/local/edit_gtf/main"
include {R as ARCHR_PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_ArchR_preprocessing.R", checkIfExists: true) )
include {R as ARCHR_FILTERING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/2_ArchR_filtering.R", checkIfExists: true) )
include {R as ARCHR_CLUSTERING_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/3_ArchR_clustering.R", checkIfExists: true) )
include {R as ARCHR_CLUSTERING_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/3_ArchR_clustering.R", checkIfExists: true) )
include {R as ARCHR_FILTER_CLUSTERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/4_ArchR_filter_clusters.R", checkIfExists: true) )

workflow ARCHR_PROCESSING {
    take:
    input

    main:
    input // [[meta], [cellranger_output, galgal_reference]]
        .set {ch_input}

    EDIT_GTF ( input ) //edits the gtf file to add 'chr' to chromosome names

    EDIT_GTF.out 
        .combine(ch_input) //[[meta], temp.gtf, [meta], cellranger_output, galgal_reference]]
        .map{[it[0], it[[1]] + it[3]]} //[[meta], [temp.gtf, cellranger_output, galgal_reference]]
        .set {ch_input_modified} // ch_metadata: [[meta], [cellranger_output, gtf]]

    ARCHR_PREPROCESSING( ch_input_modified )
    // creates arrow files and ArchR project filtered with generous thresholds
    ARCHR_FILTERING( ARCHR_PREPROCESSING.out )
    // filters whole data globally (+ on a per sample basis?)
    ARCHR_CLUSTERING_PREFILTER( ARCHR_FILTERING.out )
    // iterative LSI, cluster - non deterministic!
    ARCHR_FILTER_CLUSTERS( ARCHR_CLUSTERING_PREFILTER.out )
    // filters poor quality clusters from whole dataset
    ARCHR_CLUSTERING_POSTFILTER( ARCHR_FILTER_CLUSTERS.out )
    // iterative LSI, cluster - non deterministic!

    // add a script here that makes some plots with gene scores for full dataset

    //emit full filtered and clustered dataset:
    archr_filtered_full = ARCHR_FILTER_CLUSTERS.out
}
