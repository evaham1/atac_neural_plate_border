#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as ARCHR_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_integration.R", checkIfExists: true) )

workflow ARCHR_INTEGRATION {
    take:
    input_ch

    main:
    ARCHR_INTEGRATE ( input_ch )
    
    // script to integrate datasets, run on full data and on subsets


    //emit integrated ArchR objects:
    emit:
    archr_integrated_full = ARCHR_INTEGRATE.out
}
