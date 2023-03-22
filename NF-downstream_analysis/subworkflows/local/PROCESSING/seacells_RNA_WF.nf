#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

///// Scripts to run SEACells computation and re-run processing on metacells
// Convert to Anndata
include {R as SEURAT_TO_ANNDATA} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )
// Run SEACells
include {PYTHON as CALCULATE_SEACELLS} from "$moduleDir/../../../modules/local/python/main"               addParams(script: file("$moduleDir/../../../bin/seacells/SEACells_computation.py", checkIfExists: true) )
// Re-process SEACells in R
include {R as META_TO_SEURAT} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seacells_meta_to_seurat.R", checkIfExists: true) )
include {R as CLASSIFY_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/RNA_processes/state_classification.R", checkIfExists: true) )
// Convert back to Anndata
include {R as SEURAT_TO_ANNDATA_PROCESSED} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow SEACELLS_RNA_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]
    ch_BNM

    main:

    input.set { ch_seurat }

    // Convert seurat to Anndata object
    SEURAT_TO_ANNDATA( input )

    //////// Run SEACells /////////
    CALCULATE_SEACELLS( SEURAT_TO_ANNDATA.out ) // Python script to calculate seacells on AnnData object

    CALCULATE_SEACELLS.out.view()

    // Process resulting metacells - need to input original seurat object and the anndata exported data
    ch_combined = CALCULATE_SEACELLS.out
            .combine(ch_seurat)
            .map{[it[0], it[1] + it[3]]}
    //ch_combined.view() //[[sample_id:Test], [./plots, ./rds_files, ./ss8_clustered_data.RDS]]
    META_TO_SEURAT( ch_combined ) // this outputs 2 seurat objects, one full object with metacell assignments added and one summarised seurat
    
    // Filter output so only take summarised seurat object??

    // Re-run cell state classification on metacells
    ch_state_classification = META_TO_SEURAT.out
        .combine(ch_BNM)
        .map{[it[0], it[1] + it[3]]}

    ch_state_classification.view()

    CLASSIFY_METACELLS( ch_state_classification )

    // Convert back into Anndata object for integration
    //SEURAT_TO_ANNDATA_PROCESSED( ch )

    emit:
    anndata = CALCULATE_SEACELLS.out
    seurat = META_TO_SEURAT.out
    //processed_anndata = SEURAT_TO_ANNDATA_PROCESSED.out

}
