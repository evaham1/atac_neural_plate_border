#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_gene_scores.R", checkIfExists: true) )
include {R as HEATMAP_GEX} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow CLUSTERING {
    take:
    input

    main:
    // Cluster data
    CLUSTER( input )

    // Calculate gene scores and plot heatmaps of marker genes
    GENE_SCORES( CLUSTER.out )
    HEATMAP_GEX( CLUSTER.out )

    //emit clustered and peak-called objects:
    emit:
    output = GENE_SCORES.out
}