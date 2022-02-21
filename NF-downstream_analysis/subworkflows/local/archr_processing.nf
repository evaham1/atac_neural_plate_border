#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as ARCHR_PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_preprocessing_ArchR.R", checkIfExists: true) )

workflow PROCESSING {
    take:
    input

    main:
    input
        .set { ch_fragments }

    ARCHR_PREPROCESSING( input )


    //emit:
    //signac_predicted_gex = GEX_FILTERING.out
}
