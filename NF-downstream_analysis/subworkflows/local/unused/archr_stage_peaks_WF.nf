#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// processing modules
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_peak_calling.R", checkIfExists: true) )
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow ARCHR_STAGE_PEAKS_WF {
    take:
    input

    main:

    // Re-cluster individual stages so get all the extra plots
    CLUSTER( input )

    // Call peaks on individual stages
    PEAK_CALL( CLUSTER.out )

    //outputs:
    emit:
    output = PEAK_CALL.out
}
