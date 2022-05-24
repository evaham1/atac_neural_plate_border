#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_full_data.R", checkIfExists: true) )

include {R as CLUSTER_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

workflow FULL_PROCESSING {
    take:
    input

    main:

    // filter using cell ids from stages
    FILTER_FULL( input )

    /// processing ///
    CLUSTER_POSTFILTER( FILTER_FULL.out )

    //emit filtered and clustered stage objects:
    emit:
    output = CLUSTER_POSTFILTER.out
}
