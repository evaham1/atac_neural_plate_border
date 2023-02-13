#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R scripts to filter peaks and then cluster them using Antler
include {R as FILTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/cluster_peaks/1_filter_peaks.R", checkIfExists: true) )
include {R as CLUSTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/cluster_peaks/2_antler_calculate_peak_modules.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// cluster peaks

workflow CLUSTER_PEAKS_WF {
    take:
    input //should be TransferLabels object with metacells labelled + summarise_counts.csv

    main:

    input.view()
    
    // Filter peaks based on annotation and variability
    FILTER_PEAKS( input )

    // Cluster peaks using Antler package
    CLUSTER_PEAKS( FILTER_PEAKS.out )

    emit:
    test_output = input
}
