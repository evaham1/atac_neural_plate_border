#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R and Python scripts to run SEACells computation
include {R as EXPORT_DATA_FOR_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/1_export_data_from_ArchR.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// clustering peaks together

workflow CLUSTER_PEAKS {
    take:
    input_ch //should just be TransferLabels object

    main:
    
    //////// Run SEACells /////////
    EXPORT_DATA_FOR_SEACELLS( input_ch ) // R script to export data to run seacells computation


    emit:
    test_output = EXPORT_DATA_FOR_SEACELLS.out
}
