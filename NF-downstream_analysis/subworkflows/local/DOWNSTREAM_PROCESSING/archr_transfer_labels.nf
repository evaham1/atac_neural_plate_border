#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//////////////////        Transfer labels       ///////////////////
// peak call and visualise on full dataset (transfer labels object)
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/ArchR_utilities/transfer_labels.R", checkIfExists: true) )
include {R as PEAK_CALL_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as HEATMAP_PEAKS_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )
include {R as HEATMAP_GEX_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )

// explore peak changes over time
include {R as DIFF_PEAKS_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/diff_peaks_stages.R", checkIfExists: true) )
include {R as DIFF_PEAKS_CLUSTERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/diff_peaks_clusters.R", checkIfExists: true) )

// cluster peaks into modules
include {R as CLUSTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/Clustering_peaks.R", checkIfExists: true) )

//////////////////        Look for enhancers       ///////////////////
// look for enhancers
// include {R as SE_CALCULATE_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/calculate_se.R", checkIfExists: true) )
// include {R as FINDING_ENHANCERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Finding_enhancers/finding_enhancers.R", checkIfExists: true) )

// // plots for enhancers
// include {R as PLOT_MANUALLY_FILTERED_ENHANCERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_manually_filtered_enhancers.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow TRANSFER_LABELS {
    take:
    input_ch    // [[sample_id:HH5], [/rds_files/HH5_Save-ArchR]]
                // [[sample_id:HH6], [/rds_files/HH6_Save-ArchR]]
                // etc

    main:

    ch_combined = input_ch
        .map{meta, output -> output}
        .collect()
        .map { output -> [[sample_id:'FullData'], output] } 
        .view() // [[sample_id:FullData], [/rds_files/HH5_Save-ArchR, /rds_files/HH6_Save-ArchR, /rds_files/HH7_Save-ArchR, /rds_files/ss4_Save-ArchR, /rds_files/ss8_Save-ArchR, /rds_files/FullData_Save-ArchR] ]


    //debugging
    ch_combined_test = ch_combined.map{meta, output -> [meta, output[0]]}
    ch_combined_test.view() //[[sample_id:FullData], /rds_files/HH5_Save-ArchR]
    TRANSFER_LABELS( ch_combined_test )
    
    // visualise clusters from individual stages on full dataset
    // THIS IS WHERE THE STACK OVERFLOW ERROR IS:
    //TRANSFER_LABELS( ch_combined ) // transfers cluster labels from stage data onto full data

    // HEATMAP_GEX_TL( TRANSFER_LABELS.out )
    // PEAK_CALL_TL( TRANSFER_LABELS.out )
    // HEATMAP_PEAKS_TL( PEAK_CALL_TL.out )

    // // visualise differential peaks across full data
    // DIFF_PEAKS_STAGES( PEAK_CALL_TL.out )
    // DIFF_PEAKS_CLUSTERS( PEAK_CALL_TL.out )

    // // cluster peaks into modules
    // CLUSTER_PEAKS( PEAK_CALL_TL.out )


    ////////////////////////////////////////////////////////////////////////////////////////////
                    /// FINDING ENHANCERS   ///

    // finding enhancers
//     SE_CALCULATE_TL( PEAK_CALL_TL.out )
//     FINDING_ENHANCERS( SE_CALCULATE_TL.out )
//     PLOT_MANUALLY_FILTERED_ENHANCERS( SE_CALCULATE_TL.out )

    emit:
    //transfer_label_peaks = PEAK_CALL_TL.out
    ch_combined_output = ch_combined

}
