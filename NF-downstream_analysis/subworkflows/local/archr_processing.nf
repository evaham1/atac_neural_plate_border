#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as ARCHR_PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_ArchR_preprocessing.R", checkIfExists: true) )
include {R as ARCHR_FILTERING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_ArchR_filtering.R", checkIfExists: true) )

workflow ARCHR_PROCESSING {
    take:
    input

    main:
    input
        .set { ch_fragments }

    ARCHR_PREPROCESSING( input )
    ARCHR_FILTERING( ARCHR_PREPROCESSING.out )

    //emit:
    //signac_predicted_gex = GEX_FILTERING.out
}
