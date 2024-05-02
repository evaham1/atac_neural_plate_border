#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//////////////////        Transfer labels       ///////////////////
// transfer labels
include {R as TRANSFER_LABELS_NEW} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/ArchR_utilities/transfer_labels.R", checkIfExists: true) )

// cluster
include {R as CLUSTER_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )

// call peaks
include {R as PEAK_CALL_TL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow TRANSFER_LABELS {
    take:
    input_ch

    main:

    // transfer over the cluster and integrated labels from stages onto filtered full data
    // for cell ids in full data which are not in stage data (ie contamination) - filter them
    TRANSFER_LABELS_NEW( input_ch )

    // recluster data now that some cells have been removed
    CLUSTER_TL( TRANSFER_LABELS_NEW.out )

    // call peaks
    PEAK_CALL_TL( CLUSTER_TL.out )

    emit:
    transfer_label_peaks = PEAK_CALL_TL.out

}
