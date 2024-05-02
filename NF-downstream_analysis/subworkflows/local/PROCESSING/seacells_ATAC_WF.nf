#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Convert to Anndata
include {R as ARCHR_EXPORT_DATA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/Export_data_from_ArchR.R", checkIfExists: true) )
include {PYTHON as CREATE_ANNDATA} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/data_conversion/ATAC_exports_to_AnnData.py", checkIfExists: true) )
// Run SEACells
include {PYTHON as CALCULATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/seacells/SEACells_computation.py", checkIfExists: true) )
// Re-process SEACells in R
include {R as META_TO_SEURAT_ATAC} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seacells_meta_to_seurat_ATAC.R", checkIfExists: true) )
include {R as PROCESS_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/seacells/process_seurat_ATAC.R", checkIfExists: true) )
include {R as CLASSIFY_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/seacells/state_classification.R", checkIfExists: true) )
// Convert back to Anndata
include {R as SEURAT_TO_ANNDATA_PROCESSED_ATAC} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )
// Rename SEACell outputs
include {R as RENAME_SEACELL_OUTPUTS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/SEACell_exports_rename.R", checkIfExists: true) )

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

    //ARCHR_EXPORT_DATA.out.view()
// [[sample_id:HH6], [ArchRLogs, exported_ArchR_data]]
// [[sample_id:HH5], [ArchRLogs, exported_ArchR_data]]
// [[sample_id:HH7], [ArchRLogs, exported_ArchR_data]]
// [[sample_id:ss4], [ArchRLogs, exported_ArchR_data]]
// [[sample_id:ss8], [ArchRLogs, exported_ArchR_data]]
    
    //////// Run SEACells /////////
    CALCULATE_SEACELLS( CREATE_ANNDATA.out ) // Python script to calculate seacells on AnnData object

    //CALCULATE_SEACELLS.out.view()
// [[sample_id:HH6], [exported_data, plots, rds_files]]
// [[sample_id:HH5], [exported_data, plots, rds_files]]
// [[sample_id:HH7], [exported_data, plots, rds_files]]
// [[sample_id:ss4], [exported_data, plots, rds_files]]
// [[sample_id:ss8], [exported_data, plots, rds_files]]

    //////// Convert to Seurat using gene scores matrix /////////
    CALCULATE_SEACELLS.out
            .concat( ARCHR_EXPORT_DATA.out )
            .groupTuple( by:0 )
            //.view() //[[sample_id:HH6], [[exported_data, plots, rds_files], [ArchRLogs, exported_ArchR_data]]]
            .map{ meta, data -> [meta, [data[0][0], data[1][1]]]}
            .set {ch_combined}
    //ch_combined.view()
// [[sample_id:HH5], [exported_data, exported_ArchR_data]]
// [[sample_id:HH6], [exported_data, exported_ArchR_data]]
// [[sample_id:HH7], [exported_data, exported_ArchR_data]]
// [[sample_id:ss4], [exported_data, exported_ArchR_data]]
// [[sample_id:ss8], [exported_data, exported_ArchR_data]]
    META_TO_SEURAT_ATAC( ch_combined ) // input needs to be ArchR exported gene score matrix (in 'exported_ArchR_data') and cell_metadata.csv (in 'exported_data') from seacells computation

    //////// Process metacells Seurat object /////////
    PROCESS_METACELLS( META_TO_SEURAT_ATAC.out )

    // // Run cell state classification on metacells
    // ch_state_classification = PROCESS_METACELLS.out
    //     .combine(ch_BNM)
    //     .map{[it[0], it[1] + it[2]]}
    // //ch_state_classification.view() //[[sample_id:ss4], [plots, rds_files, binary_knowledge_matrix_contam.csv]]
    // CLASSIFY_METACELLS( ch_state_classification )

    //////// Convert to Anndata /////////
    SEURAT_TO_ANNDATA_PROCESSED_ATAC( PROCESS_METACELLS.out )

    //////// Rename SEACell outputs for downstream peak modules /////////
    RENAME_SEACELL_OUTPUTS( CALCULATE_SEACELLS.out )

    emit:
    seacells_anndata = CALCULATE_SEACELLS.out
    seacells_seurat_objects = META_TO_SEURAT_ATAC.out
    seacells_seurat_processed = PROCESS_METACELLS.out
    //seacells_seurat_processed_classified = CLASSIFY_METACELLS.out
    seacells_anndata_processed_classified = SEURAT_TO_ANNDATA_PROCESSED_ATAC.out
    seacell_outputs_named = RENAME_SEACELL_OUTPUTS.out

}
