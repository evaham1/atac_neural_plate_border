#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PREPROCESS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_preprocessing.R", checkIfExists: true) )

workflow PREPROCESSING {
    take:
    input

    main:
    PREPROCESS( input )

    //emit full filtered and clustered dataset:
    emit:
    output = PREPROCESS.out
}