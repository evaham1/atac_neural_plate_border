#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

///// Scripts to run SEACells computation and re-run processing on metacells
// Convert to Anndata
include {R as SEURAT_TO_ANNDATA} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )
// Run SEACells
include {PYTHON as CALCULATE_SEACELLS} from "$moduleDir/../../../modules/local/python/main"               addParams(script: file("$moduleDir/../../../bin/seacells/SEACells_computation.py", checkIfExists: true) )
// Re-process SEACells in R
include {R as META_TO_SEURAT_RNA} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seacells_meta_to_seurat_RNA.R", checkIfExists: true) )
include {R as PROCESS_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/seacells/process_seurat_RNA.R", checkIfExists: true) )
include {R as CLASSIFY_METACELLS} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/seacells/state_classification.R", checkIfExists: true) )
// Convert back to Anndata
include {R as SEURAT_TO_ANNDATA_PROCESSED_RNA} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )
// Transfer RNA metacell labels onto full data object to check how well labels match
include {R as TRANSFER_METACELL_LABELS_RNA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/Transfer_metacell_labels_RNA.R", checkIfExists: true) )

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

    //////// Convert Seurat to AnnData /////////
    SEURAT_TO_ANNDATA( input )

    //////// Run SEACells /////////
    CALCULATE_SEACELLS( SEURAT_TO_ANNDATA.out ) // Python script to calculate seacells on AnnData object

    //CALCULATE_SEACELLS.out.view()
    // [[sample_id:HH6], [./exported_data, ./plots, ./rds_files]]

    //////// Convert AnnData to Seurat /////////
    CALCULATE_SEACELLS.out
            .concat( ch_seurat )
            .groupTuple( by:0 )
            .map{ meta, data -> [meta, [data[0][0], data[1][0]]]}
            .set { ch_combined }

    META_TO_SEURAT_RNA( ch_combined ) // this outputs 2 seurat objects, one full object with metacell assignments added and one summarised seurat

    // extract just the full 
    META_TO_SEURAT_RNA.out
        .map{ [it[0], [it[1].findAll{it =~ /rds_files_full/}[0].listFiles()[0]]] }
        .set{ metacell_assignments }
    metacell_assignments.view()
    //     [[sample_id:HH5], seacells_seurat.RDS]
    //     [[sample_id:ss4], seacells_seurat.RDS]
    //     [[sample_id:HH6], seacells_seurat.RDS]
    //     [[sample_id:ss8], seacells_seurat.RDS]
    //     [[sample_id:HH7], seacells_seurat.RDS]
                // [[sample_id:HH5], [[seacells_seurat.RDS, seurat.RDS]]]
                // [[sample_id:HH6], [[seacells_seurat.RDS, seurat.RDS]]]
                // [[sample_id:ss4], [[seacells_seurat.RDS, seurat.RDS]]]
                // [[sample_id:ss8], [[seacells_seurat.RDS, seurat.RDS]]]
                // [[sample_id:HH7], [[seacells_seurat.RDS, seurat.RDS]]]

    //////// Process metacells Seurat object /////////
    PROCESS_METACELLS( META_TO_SEURAT_RNA.out )

    // Re-run cell state classification on metacells
    ch_state_classification = PROCESS_METACELLS.out
        .combine(ch_BNM)
        .map{[it[0], it[1] + it[2]]}
    //ch_state_classification.view()
    CLASSIFY_METACELLS( ch_state_classification )

    //////// Convert to Anndata /////////
    SEURAT_TO_ANNDATA_PROCESSED_RNA( CLASSIFY_METACELLS.out )

    /// Transfer metacell labels onto single cells and check how they compare ///
    CLASSIFY_METACELLS.out
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( metacell_assignments )
            .groupTuple( by:0 )
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            .set { ch_combined_2 }
    TRANSFER_METACELL_LABELS_RNA( ch_combined_2 )
    


    emit:
    seacells_anndata = CALCULATE_SEACELLS.out
    seacells_seurat_objects = META_TO_SEURAT_RNA.out
    seacells_seurat_processed = PROCESS_METACELLS.out
    seacells_seurat_processed_classified = CLASSIFY_METACELLS.out
    seacells_anndata_processed_classified = SEURAT_TO_ANNDATA_PROCESSED_RNA.out

}
