#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Python scripts to run SEACells computation
include {R as ARCHR_EXPORT_DATA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/Export_data_from_ArchR.R", checkIfExists: true) )

include {PYTHON as CREATE_ANNDATA} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/data_conversion/ATAC_exports_to_AnnData.py", checkIfExists: true) )
include {PYTHON as CALCULATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/SEACells_computation.py", checkIfExists: true) )
include {PYTHON as EXPORT_DATA_FROM_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/data_conversion/Export_data_from_AnnData.py", checkIfExists: true) )

//include {R as CHECK_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/Check_seacell_purity.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow SEACELLS_ATAC_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]

    main:

    // Convert ArchR to Anndata object
    ARCHR_EXPORT_DATA( input ) // R script to export data to run seacells computation
    CREATE_ANNDATA( ARCHR_EXPORT_DATA.out ) // Python script to read exported data into an Anndata object
    
    //////// Run SEACells /////////
    CALCULATE_SEACELLS( CREATE_ANNDATA.out ) // Python script to calculate seacells on AnnData object
    EXPORT_DATA_FROM_SEACELLS( CALCULATE_SEACELLS.out ) //Python script to export data from Anndata objects as .csv

    // Downstream processing of Metacells (manipulate channel to just extract summarised AnnData object for intergration)


    emit:
    anndata = CALCULATE_SEACELLS.out
    exports = EXPORT_DATA_FROM_SEACELLS.out

}
