#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_full_data.R", checkIfExists: true) )

include {R as CLUSTER_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as PEAK_CALL_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as PEAK_DIFF_POSTFILTER} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/Peak_calling/ArchR_diff_peaks.R", checkIfExists: true) )


workflow FULL_PROCESSING {
    take:
    input

    main:

    // filter using cell ids from stages
    FILTER_FULL( input )

    /// processing ///
    CLUSTER_POSTFILTER( FILTER_FULL.out )
    GENE_SCORES_POSTFILTER( CLUSTER_POSTFILTER.out )
    PEAK_CALL_POSTFILTER( CLUSTER_POSTFILTER.out )
    PEAK_DIFF_POSTFILTER( PEAK_CALL_POSTFILTER.out )

    //emit filtered and clustered stage objects:
    emit:
    output = PEAK_CALL_POSTFILTER.out
    gex = GENE_SCORES_POSTFILTER.out
    differential_peaks = PEAK_DIFF_POSTFILTER.out
}
