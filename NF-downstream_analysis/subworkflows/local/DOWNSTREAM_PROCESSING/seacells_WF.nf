#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Python scripts to run SEACells computation
include {PYTHON as CALCULATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/3_SEACells_computation.py", checkIfExists: true) )
include {PYTHON as EXPORT_DATA_FROM_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/4_export_data_from_AnnData.py", checkIfExists: true) )
include {R as CHECK_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/5_check_seacell_purity.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow SEACELLS_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]

    main:
    
    //////// Run SEACells /////////
    CALCULATE_SEACELLS( input ) // Python script to calculate seacells on AnnData object
    EXPORT_DATA_FROM_SEACELLS( CALCULATE_SEACELLS.out ) //Python script to export data from Anndata objects as .csv

    // //////// Check SEACells and add labels to TL object /////////
    // EXPORT_DATA_FROM_SEACELLS.out
    //         .combine(input)
    //         .map{[it[0], it[1] + it[3]]}
    //         //.view() //[[sample_id:TransferLabels], [plots, rds_files, TransferLabels_Save-ArchR]]
    //         .set {ch_combined} // combine the transferlabels object and the output from seacells calculation
    // CHECK_SEACELLS( ch_combined )

    emit:
    seacells_output = CALCULATE_SEACELLS.out
    //seacells_output_combined = CHECK_SEACELLS.out

}
