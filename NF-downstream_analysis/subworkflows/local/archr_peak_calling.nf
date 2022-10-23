#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as HEATMAP_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow PEAK_CALLING {
    take:
    input

    main:
    // Run peak calling and examine resulting differential peaks
    PEAK_CALL( CLUSTER.out )
    HEATMAP_PEAKS( PEAK_CALL.out )

    //emit clustered and peak-called objects:
    emit:
    output = PEAK_CALL.out
}