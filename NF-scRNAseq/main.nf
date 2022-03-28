#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/downstream
========================================================================================
    Github : https://github.com/nf-core/downstream
    Website: https://nf-co.re/downstream
    Slack  : https://nfcore.slack.com/channels/downstream
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { METADATA } from "$baseDir/subworkflows/local/metadata"

include {R as SPLIT} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/split_seurat.R", checkIfExists: true) )
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/subset_cluster.R", checkIfExists: true) )
include {R as STATE_CLASSIFICATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/state_classification_contam.R", checkIfExists: true) )

include {R as MERGE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/merge.R", checkIfExists: true) )
include {R as CLUSTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/subset_cluster.R", checkIfExists: true) )
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/transfer_labels.R", checkIfExists: true) )

include {R as INPUT_CHECK} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/check.R", checkIfExists: true) )
include {R as OUTPUT_CHECK} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/check.R", checkIfExists: true) )

//
// SET CHANNELS
//

// Set channel for binary knowledge matrix for cell state classification
Channel
    .value("$baseDir/binary_knowledge_matrix_contam.csv")
    .set{ch_binary_knowledge_matrix}


//
// WORKFLOW: Run main nf-core/downstream analysis pipeline
//
workflow NFCORE_DOWNSTREAM {

    METADATA( params.sample_sheet )
    //METADATA.out:
    //[[sample_id:NF-scRNA-input], [./cell_cycle_data.RDS]]

    // make some plots to check the input is exactly as expected (cell_cycle_data.RDS)
    INPUT_CHECK( METADATA.out )

    // split the cell_cycle_data object into individual stages
    SPLIT( METADATA.out )

    // filter out the HH4 stage and split remaining stages into individual channels
    SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .filter{!(it =~ /HH4/)}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }

    // for each stage cluster and classify cell states using BNM that includes contaminating populations
    CLUSTER( ch_split_run )
    CLUSTER.out
        .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]}
        .combine(ch_binary_knowledge_matrix) // Combine with binary knowledge matrix
        .map{ row -> [row[0], [row[1], row[2]]]}
        .set { ch_state_classification }    //Channel: [[meta], [rds_file, csv]]
    STATE_CLASSIFICATION( ch_state_classification )

    // Collect rds files from all stages and merge into one combined object
    ch_combined = STATE_CLASSIFICATION.out
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'all_stages'], it] } // [[meta], [rds1, rds2, rds3, ...]]
        .view()
    MERGE( ch_combined )

    // run clustering on merged data
    CLUSTER_FULL( MERGE.out )

    // // Transfer labels from individual stages to merged data
    // ch_labels = STATE_CLASSIFICATION.out
    //     .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
    //     .collect()
    //     .map { [[sample_id:'all_stages'], it] } // [[meta], [rds1, rds2, rds3, ...]]
    //     .combine( CLUSTER_FULL.out ) //[[sample_id:all_stages], [HH7, ss8, HH6, ss4, HH4, HH5], [sample_id:NF-scRNA-input], [rds_files, plots]]
    //     .map{[it[0], it[1] + it[3]]} //[[sample_id:all_stages], [HH7, ss8, HH6, ss4, HH4, HH5, rds_files, plots]
    //     //.view() //[[sample_id:all_stages], [HH6, HH4, ss8, ss4, HH7, HH5, cell_cycle_data.RDS]]
    // TRANSFER_LABELS( ch_labels )
}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_DOWNSTREAM ()
}


/*
========================================================================================
    THE END
========================================================================================
*/