#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R and Python scripts to convert ArchR data to Anndata object
include {R as ARCHR_EXPORT_DATA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/1_export_data_from_ArchR.R", checkIfExists: true) )
include {PYTHON as CREATE_ANNDATA} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/2_exports_to_AnnData.py", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow ARCHR_TO_ANNDATA_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]

    main:
    ARCHR_EXPORT_DATA( input ) // R script to export data to run seacells computation
    CREATE_ANNDATA( ARCHR_EXPORT_DATA.out ) // Python script to read exported data into an Anndata object

    emit:
    anndata = CREATE_ANNDATA.out
    exported_archr = ARCHR_EXPORT_DATA.out

}
