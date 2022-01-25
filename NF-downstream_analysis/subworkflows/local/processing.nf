#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_preprocessing.R", checkIfExists: true) )
include {R as FILTERING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/2_filtering.R", checkIfExists: true) )
include {R as GENE_ACTIVITY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/3_gene_activity.R", checkIfExists: true) )

workflow PROCESSING {
    take:
    input

    main:
    PREPROCESSING( input )
    FILTERING(PREPROCESSING.out)
    GENE_ACTIVITY(FILTERING.out)

    //emit:
    //preprocessing_out = PREPROCESSING.out
}
