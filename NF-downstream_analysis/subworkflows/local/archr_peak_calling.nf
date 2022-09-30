#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_gene_scores.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as HEATMAP_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )
include {R as HEATMAP_GEX} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow PEAK_CALLING {
    take:
    input

    main:
    
    CLUSTER( input )
    
    GENE_SCORES( CLUSTER.out )
    HEATMAP_GEX( CLUSTER.out )

    PEAK_CALL( CLUSTER.out )
    HEATMAP_PEAKS( PEAK_CALL.out )

    //emit clustered and peak-called objects:
    emit:
    output = PEAK_CALL.out
}
