#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//////////////////        NPB subset       ///////////////////
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

workflow NPB_SUBSET {
    take:
    input_ch

    main:


    ////////////////////////////////////////////////////////////////////////////////////////////
                    /// NPB SUBSET   ///

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
    //transfer_label_peaks = PEAK_CALL_TL.out
    //npb_peaks = PEAK_CALL_NPB.out
}
