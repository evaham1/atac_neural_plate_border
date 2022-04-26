#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER_NONE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as FILTER_MED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as FILTER_HIGH} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )

workflow QC_EXPLORATION {
    take:
    input

    main:

    /// FILTER NONE ///
    FILTER_NONE( input )
    
    /////////////////////////
    
    



    //emit full filtered and clustered dataset:
    emit:
    output = FILTER.out
}