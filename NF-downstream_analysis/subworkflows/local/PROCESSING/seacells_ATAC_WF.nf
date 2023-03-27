#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Convert to Anndata
include {R as ARCHR_EXPORT_DATA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/Export_data_from_ArchR.R", checkIfExists: true) )
include {PYTHON as CREATE_ANNDATA} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/data_conversion/ATAC_exports_to_AnnData.py", checkIfExists: true) )
// Run SEACells
include {PYTHON as CALCULATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/SEACells_computation.py", checkIfExists: true) )
// Re-process SEACells in R
include {R as META_TO_SEURAT_ATAC} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seacells_meta_to_seurat_ATAC.R", checkIfExists: true) )
include {R as PROCESS_METACELLS_ATAC} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/Metacell_processes/process_seurat_ATAC.R", checkIfExists: true) )
include {R as CLASSIFY_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/Metacell_processes/state_classification.R", checkIfExists: true) )
// Convert back to Anndata
include {R as SEURAT_TO_ANNDATA_PROCESSED} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow SEACELLS_ATAC_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]
    ch_BNM
    //RNA_metacells // processed seurat object with RNA metacells?? maybe not needed

    main:

    //////// Convert ArchR to AnnData /////////
    ARCHR_EXPORT_DATA( input ) // R script to export data to run seacells computation
    CREATE_ANNDATA( ARCHR_EXPORT_DATA.out ) // Python script to read exported data into an Anndata object
    
    //////// Run SEACells /////////
    CALCULATE_SEACELLS( CREATE_ANNDATA.out ) // Python script to calculate seacells on AnnData object

    //////// Convert to Seurat using gene scores matrix /////////
    ch_combined = ARCHR_EXPORT_DATA.out
        .combine(CALCULATE_SEACELLS.out)
        .map{[it[0], it[1] + it[3]]}
    ch_combined.view()
    META_TO_SEURAT_ATAC( ch_combined ) // input needs to be ArchR exported gene score matrix and cell_metadata.csv from seacells computation

    //////// Process metacells Seurat object /////////
    PROCESS_METACELLS_ATAC( META_TO_SEURAT_ATAC.out )

    // Run cell state classification on metacells
    ch_state_classification = PROCESS_METACELLS_ATAC.out
        .combine(ch_BNM)
        .map{[it[0], it[1] + it[2]]}
    ch_state_classification.view()
    CLASSIFY_METACELLS( ch_state_classification )

    //////// Convert to Anndata /////////
    SEURAT_TO_ANNDATA_PROCESSED( CLASSIFY_METACELLS.out )

    emit:
    seacells_anndata = CALCULATE_SEACELLS.out
    seacells_seurat_objects = META_TO_SEURAT_ATAC.out
    seacells_seurat_processed = PROCESS_METACELLS_ATAC.out
    seacells_seurat_processed_classified = CLASSIFY_METACELLS.out
    seacells_anndata_processed_classified = SEURAT_TO_ANNDATA_PROCESSED.out

}
