#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R and Python scripts to run SEACells computation
include {R as EXPORT_DATA_FOR_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/1_export_data_from_ArchR.R", checkIfExists: true) )
include {PYTHON as CREATE_ANNDATA} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/2_exports_to_AnnData.py", checkIfExists: true) )
include {PYTHON as CALCULATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/3_SEACells_computation.py", checkIfExists: true) )
include {PYTHON as EXPORT_DATA_FROM_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/4_export_data_from_AnnData.py", checkIfExists: true) )
include {R as CHECK_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/5_check_seacell_purity.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow CALCULATE_SEACELLS_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]

    main:
    
    //////// Run SEACells /////////
    EXPORT_DATA_FOR_SEACELLS( input ) // R script to export data to run seacells computation
    CREATE_ANNDATA( EXPORT_DATA_FOR_SEACELLS.out ) // Python script to read exported data into an Anndata object
    CALCULATE_SEACELLS( CREATE_ANNDATA.out ) // Python script to calculate seacells on AnnData object
    EXPORT_DATA_FROM_SEACELLS( CALCULATE_SEACELLS.out ) //Python script to export data from Anndata objects as .csv

    //////// Check SEACells and add labels to TL object /////////
    EXPORT_DATA_FROM_SEACELLS.out
            .combine(input_ch)
            .map{[it[0], it[1] + it[2]]}
            .view()
            .set {ch_combined} // combine the transferlabels object and the output from seacells calculation
    CHECK_SEACELLS( ch_combined )

    emit:
    seacells_output = CALCULATE_SEACELLS.out
    seacells_output_combined = CHECK_SEACELLS.out

}
