#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R and Python scripts to run SEACells computation
include {R as EXPORT_DATA_FOR_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/1_export_data_from_ArchR.R", checkIfExists: true) )
include {PYTHON as CREATE_ANNDATA} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/2_exports_to_AnnData.py", checkIfExists: true) )
include {PYTHON as CALCULATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/2_exports_to_AnnData.py", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// clustering peaks together

workflow CLUSTER_PEAKS {
    take:
    input_ch //should just be TransferLabels object

    main:
    
    //////// Run SEACells /////////
    EXPORT_DATA_FOR_SEACELLS( input_ch ) // R script to export data to run seacells computation
    CREATE_ANNDATA( EXPORT_DATA_FOR_SEACELLS.out ) // Python script to read exported data into an Anndata object
    //CALCULATE_SEACELLS( CREATE_ANNDATA ) // Python script to calculate seacells on AnnData object

    emit:
    test_output = CREATE_ANNDATA.out
}
