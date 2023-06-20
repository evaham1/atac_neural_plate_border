#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// integration
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_peak_calling.R", checkIfExists: true) )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow ARCHR_STAGE_PEAKS_WF {
    take:
    input

    main:

    // Re-cluster individual stages so get all the extra plots
    CLUSTER( input )

    // Call peaks on individual stages
    PEAK_CALL ( CLUSTER.out )

    // Find differential peaks between clusters and see how these change over time

    //outputs:
    emit:
    peak_called = PEAK_CALL.out
}
