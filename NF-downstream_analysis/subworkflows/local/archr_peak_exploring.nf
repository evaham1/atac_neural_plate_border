#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// compare stages
//include {R as COMPARE_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/compare_stages.R", checkIfExists: true) )

// visualise on full dataset (transfer labels object)
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/ArchR_utilities/transfer_labels.R", checkIfExists: true) )
include {R as PEAK_CALL_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as HEATMAP_PEAKS_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )
include {R as HEATMAP_GEX_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Visualisations/plot_marker_heatmaps.R", checkIfExists: true) )

// explore peak changes over time
include {R as DIFF_PEAKS_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/diff_peaks_stages.R", checkIfExists: true) )
include {R as DIFF_PEAKS_CLUSTERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/diff_peaks_clusters.R", checkIfExists: true) )

// look for enhancers
include {R as SE_CALCULATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/calculate_se.R", checkIfExists: true) )
include {R as FINDING_ENHANCERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/finding_enhancers.R", checkIfExists: true) )

// // look for enhancers in just NPB cells
// include {R as SUBSET_NPB} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
// include {R as CLUSTER_NPB} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
// include {R as PEAK_CALL_NPB} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )

// // just NPB at each stage
// include {R as SUBSET_NPB_HH7} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
// include {R as CLUSTER_NPB_HH7} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

// include {R as SUBSET_NPB_SS4} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
// include {R as CLUSTER_NPB_SS4} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

// include {R as SUBSET_NPB_SS8} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
// include {R as CLUSTER_NPB_SS8} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

// // look for enhancers in just NPB cells - subset using just labels
// include {R as SUBSET_NPB_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
// include {R as CLUSTER_NPB_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
// include {R as PEAK_CALL_NPB_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow PEAK_EXPLORING {
    take:
    input_ch

    main:
    // collect all integrated rds objects into a single element in channel
    ch_combined = input_ch // Collect integrated atac objects
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            .collect()
            .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]
            .view()

    // compare variability/how many differential peaks we have at different stages
    //COMPARE_STAGES( ch_combined )

    // visualise clusters from individual stages on full dataset
    TRANSFER_LABELS( ch_combined ) // transfers cluster labels from stage data onto full data
    HEATMAP_GEX_TL( TRANSFER_LABELS.out )
    PEAK_CALL_TL( TRANSFER_LABELS.out )
    HEATMAP_PEAKS_TL( PEAK_CALL_TL.out )

    // visualise differential peaks across full data
    //DIFF_PEAKS_STAGES( PEAK_CALL_TL.out )
    //DIFF_PEAKS_CLUSTERS( PEAK_CALL_TL.out )

    // finding enhancers
    SE_CALCULATE( PEAK_CALL_TL.out )
    FINDING_ENHANCERS( SE_CALCULATE.out )

    // subset NPB
    // SUBSET_NPB( TRANSFER_LABELS.out )
    // CLUSTER_NPB( SUBSET_NPB.out )
    // PEAK_CALL_NPB( CLUSTER_NPB.out )

    // // subset NPB at each stage
    // SUBSET_NPB_HH7( TRANSFER_LABELS.out )
    // CLUSTER_NPB_HH7( SUBSET_NPB_HH7.out )
    // SUBSET_NPB_SS4( TRANSFER_LABELS.out )
    // CLUSTER_NPB_SS4( SUBSET_NPB_SS4.out )
    // SUBSET_NPB_SS8( TRANSFER_LABELS.out )
    // CLUSTER_NPB_SS8( SUBSET_NPB_SS8.out )

    // // subset NPB using just labels
    // SUBSET_NPB_LABELS( TRANSFER_LABELS.out )
    // CLUSTER_NPB_LABELS( SUBSET_NPB_LABELS.out )
    // PEAK_CALL_NPB_LABELS( CLUSTER_NPB_LABELS.out )
    
    emit:
    transfer_label_peaks = PEAK_CALL_TL.out
    //npb_peaks = PEAK_CALL_NPB.out
}
