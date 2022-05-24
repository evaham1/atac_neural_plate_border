#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// integration
include {R as UNCON_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_unconstrained_integration.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )

// remove contamination
include {R as SUBSET_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
include {R as CLUSTER_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY_FILTERED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )

// characterise gene scores
include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as HEATMAP_GEX} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/plot_marker_heatmaps.R", checkIfExists: true) )

// characterise peaks
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as PEAK_DIFF} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/Peak_calling/ArchR_diff_peaks_analysis.R", checkIfExists: true) )
include {R as HEATMAP_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/plot_marker_heatmaps.R", checkIfExists: true) )

// visualise full dataset
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/ArchR_preprocessing/transfer_labels.R", checkIfExists: true) )
include {R as PEAK_CALL_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as PEAK_DIFF_TL} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/Peak_calling/ArchR_diff_peaks.R", checkIfExists: true) )


workflow INTEGRATING {
    take:
    input_ch

    main:

    // Integrate full data and split stage data
    UNCON_INTEGRATE ( input_ch )

    // Label clusters based on most frequent label within each cluster
    CLUSTER_IDENTIFY ( UNCON_INTEGRATE.out )
    
    // // Filter contaminating cells from all channels and re-cluster all channels
    // SUBSET_INTEGRATION ( UNCON_INTEGRATE.out )
    // CLUSTER_INTEGRATION ( SUBSET_INTEGRATION.out )
    // CLUSTER_IDENTIFY_FILTERED ( CLUSTER_INTEGRATION.out )

    // Characterise clusters using gene scores
    GENE_SCORES( UNCON_INTEGRATE.out )
    HEATMAP_GEX( UNCON_INTEGRATE.out )

    // Characterise clusters using peaks
    PEAK_CALL( UNCON_INTEGRATE.out )
    //PEAK_DIFF( PEAK_CALL.out )
    HEATMAP_PEAKS( PEAK_CALL.out )

    // Look at the stage clusters on the full dataset
    ch_combined = UNCON_INTEGRATE.out // Collect integrated atac objects
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            .collect()
            .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]
            .view()
    TRANSFER_LABELS( ch_combined )
    PEAK_CALL_TL( TRANSFER_LABELS.out )
    PEAK_DIFF_TL( PEAK_CALL_TL.out )


    //emit integrated ArchR objects:
    emit:
    //archr_integrated_full = CLUSTER_INTEGRATION.out
    transfer_labels = TRANSFER_LABELS.out
}
