#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_preprocessing.R", checkIfExists: true) )


workflow PROCESSING {
    take:
    input

    main:
    PREPROCESSING( input )

    //emit:
    //preprocessing_out = PREPROCESSING.out
}
