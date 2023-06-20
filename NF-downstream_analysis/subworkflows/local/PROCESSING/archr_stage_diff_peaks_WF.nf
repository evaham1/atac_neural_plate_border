#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// processing modules
include {R as CALCULATE_SE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_peaks/calculate_se.R", checkIfExists: true) )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow ARCHR_STAGE_DIFF_PEAKS_WF {
    take:
    input

    main:

    // Calculate se object for each stage
    CALCULATE_SE( input )

    // Make plots of differential peaks from different annotations
    
    // Find differential peaks between clusters and see how these change over time

    // Idetify which clusters become different first?

    //outputs:
    emit:
    peak_called = PEAK_CALL.out
}
