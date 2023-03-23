#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

///// Scripts to run SEACells computation and re-run processing on metacells
// Convert to Anndata
include {R as SEURAT_TO_ANNDATA} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )
// Run SEACells
include {PYTHON as CALCULATE_SEACELLS} from "$moduleDir/../../../modules/local/python/main"               addParams(script: file("$moduleDir/../../../bin/seacells/SEACells_computation.py", checkIfExists: true) )
// Re-process SEACells in R
include {R as META_TO_SEURAT} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seacells_meta_to_seurat.R", checkIfExists: true) )
include {R as PROCESS_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/RNA_processes/process_seurat.R", checkIfExists: true) )
include {R as CLASSIFY_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/RNA_processes/state_classification.R", checkIfExists: true) )
// Convert back to Anndata
include {R as SEURAT_TO_ANNDATA_PROCESSED} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow SEACELLS_RNA_WF {
    take:
    input
    ch_BNM

    main:

    //input.view()
    // [[sample_id:HH5], [/ready_for_integration/HH5_splitstage_data/rds_files/HH5_clustered_data.RDS]]
    // [[sample_id:HH6], [/ready_for_integration/HH6_splitstage_data/rds_files/HH6_clustered_data.RDS]]
    // [[sample_id:HH7], [/ready_for_integration/HH7_splitstage_data/rds_files/HH7_clustered_data.RDS]]
    // [[sample_id:ss4], [/ready_for_integration/ss4_splitstage_data/rds_files/ss4_clustered_data.RDS]]
    // [[sample_id:ss8], [/ready_for_integration/ss8_splitstage_data/rds_files/ss8_clustered_data.RDS]]
    //ch_BNM.view()
    /// /atac_neural_plate_border/NF-downstream_analysis/binary_knowledge_matrix_contam.csv

    input.set{ch_seurat}

    // Convert seurat to Anndata object
    SEURAT_TO_ANNDATA( input )

    //////// Run SEACells /////////
    CALCULATE_SEACELLS( SEURAT_TO_ANNDATA.out ) // Python script to calculate seacells on AnnData object

    //CALCULATE_SEACELLS.out.view()
    // [[sample_id:HH6], [./exported_data, ./plots, ./rds_files]]

    // Process resulting metacells - need to input original seurat object and the anndata exported data
    
    // ch_combined = CALCULATE_SEACELLS.out
    //         .combine(ch_seurat)
    //         .map{[it[0], it[1] + it[3]]}
    // ch_combined.view() //[[sample_id:Test], [./plots, ./rds_files, ./ss8_clustered_data.RDS]]
    
    CALCULATE_SEACELLS.out
            .concat( ch_seurat )
            .groupTuple( by:0 )
            .map{ meta, data -> [meta, [data[0][0], data[1][0]]]}
            .set {ch_combined}

    ch_combined.view()

    META_TO_SEURAT( ch_combined ) // this outputs 2 seurat objects, one full object with metacell assignments added and one summarised seurat

    // Re-process summarised seurat object
    PROCESS_METACELLS( META_TO_SEURAT.out )

    // Re-run cell state classification on metacells
    ch_state_classification = PROCESS_METACELLS.out
        .combine(ch_BNM)
        .map{[it[0], it[1] + it[3]]}
    ch_state_classification.view()
    CLASSIFY_METACELLS( ch_state_classification )

    // Convert back into Anndata object for integration
    SEURAT_TO_ANNDATA_PROCESSED( CLASSIFY_METACELLS.out )

    emit:
    seacells_anndata = CALCULATE_SEACELLS.out
    seacells_seurat_objects = META_TO_SEURAT.out
    seacells_seurat_processed = PROCESS_METACELLS.out
    seacells_seurat_processed_classified = CLASSIFY_METACELLS.out
    seacells_anndata_processed_classified = SEURAT_TO_ANNDATA_PROCESSED.out

}
