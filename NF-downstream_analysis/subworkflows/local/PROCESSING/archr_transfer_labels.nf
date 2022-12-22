#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//////////////////        Transfer labels       ///////////////////
// transfer labels
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/ArchR_utilities/transfer_labels.R", checkIfExists: true) )
include {R as TRANSFER_LABELS_NEW} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/ArchR_utilities/transfer_labels_new.R", checkIfExists: true) )

// cluster
include { CLUSTERING as CLUSTERING_TL } from "$baseDir/subworkflows/local/PROCESSING/archr_clustering_and_gex"

// peak call
include {R as PEAK_CALL_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )

// visualise
include {R as HEATMAP_PEAKS_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )
include {R as HEATMAP_GEX_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )

///

// explore peak changes over time
include {R as DIFF_PEAKS_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/diff_peaks_stages.R", checkIfExists: true) )
include {R as DIFF_PEAKS_CLUSTERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/diff_peaks_clusters.R", checkIfExists: true) )

// cluster peaks into modules
include {R as CLUSTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/Clustering_peaks.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow TRANSFER_LABELS {
    take:
    input_ch

    main:

    // transfer over the cluster and integrated labels from stages onto filtered full data
    // for cell ids in full data which are not in stage data (ie contamination) - filter them - NEED TO EDIT SCRIPT TO DO THIS TOO
    TRANSFER_LABELS_NEW( input_ch )

    // recluster data now that some cells have been removed
    CLUSTERING_TL( TRANSFER_LABELS_NEW.out )

    // call peaks
    PEAK_CALL_TL( CLUSTERING_TL.out )

    // visualise the transfer_labels object 
    //HEATMAP_PEAKS_TL( PEAK_CALL_TL.out )
    //DIFF_PEAKS_STAGES( PEAK_CALL_TL.out )
    //DIFF_PEAKS_CLUSTERS( PEAK_CALL_TL.out )


    emit:
    transfer_label_peaks = PEAK_CALL_TL.out
    ch_combined_output = input_ch

}
