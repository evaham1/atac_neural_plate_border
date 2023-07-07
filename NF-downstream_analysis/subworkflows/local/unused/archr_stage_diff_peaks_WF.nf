#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// processing modules
include {R as CALCULATE_SE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_peaks/calculate_se.R", checkIfExists: true) )
include {R as PLOT_DIFF_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_peaks/diff_peaks_plots.R", checkIfExists: true) )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow ARCHR_STAGE_DIFF_PEAKS_WF {
    take:
    input

    main:

    // Calculate se object for each stage
    CALCULATE_SE( input )

    // Make plots of differential peaks for each stage (eg from different annotations, heatmaps, volcano plots, which clusters contribute to most differential peaks, etc)
    PLOT_DIFF_PEAKS( CALCULATE_SE.out )

    // Combine all stages into one channel?
    // how do these differential peaks change over time?

    //outputs:
    emit:
    output = CALCULATE_SE.out
}
