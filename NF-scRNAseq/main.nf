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
include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/transfer_labels.R", checkIfExists: true) )
include {R as SUBSET} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/subset_cells.R", checkIfExists: true) )
include {R as CLUSTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/subset_cluster.R", checkIfExists: true) )


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
    //[[sample_id:NF-scRNA-input], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/output/NF-scRNAseq/input/cell_cycle_data.RDS]]

    SPLIT( METADATA.out )

    SPLIT.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_split_run }                                                           //Channel: [[meta], rds_file]

    CLUSTER( ch_split_run )

    CLUSTER.out
        .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]}
        .combine(ch_binary_knowledge_matrix) // Combine with binary knowledge matrix
        .map{ row -> [row[0], [row[1], row[2]]]}
        .set { ch_state_classification }    //Channel: [[meta], [rds_file, csv]]

    STATE_CLASSIFICATION( ch_state_classification )

    //STATE_CLASSIFICATION output:
    //[[sample_id:HH7_splitstage_data], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/f8/2f6b0e5f0cf304fda2e33f185ee64e/plots, /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/f8/2f6b0e5f0cf304fda2e33f185ee64e/rds_files]]
    //[[sample_id:HH6_splitstage_data], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/3c/1ca35ca80a4213645f2ba0ef295ca0/plots, /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/3c/1ca35ca80a4213645f2ba0ef295ca0/rds_files]]
    //[[sample_id:HH5_splitstage_data], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/dc/2553c1c56d5fa32ed05e1634d0c373/plots, /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/dc/2553c1c56d5fa32ed05e1634d0c373/rds_files]]
    //[[sample_id:ss8_splitstage_data], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/88/4cc8853bf43cd1ad44a4898a4be9f2/plots, /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/88/4cc8853bf43cd1ad44a4898a4be9f2/rds_files]]
    //[[sample_id:HH4_splitstage_data], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/a7/201ecf787e4305b5508f95de3d6e54/plots, /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/a7/201ecf787e4305b5508f95de3d6e54/rds_files]]
    //[[sample_id:ss4_splitstage_data], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/c7/f0483ca9abe305b375f48da0cb9a46/plots, /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-scRNAseq/work/c7/f0483ca9abe305b375f48da0cb9a46/rds_files]]

    // Collect rds files from all stages
    ch_combined = STATE_CLASSIFICATION.out
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'all_stages'], it] } // [[meta], [rds1, rds2, rds3, ...]]
        .combine( METADATA.out ) //[[meta], [rds1, rds2, rds3, ...], [meta], [full.rds]]]
        .view()
        .map{[it[0], it[[1]] + it[3]]}
        .view()

    // Transfer labels from stage subsets to full data
    TRANSFER_LABELS( ch_combined )

    // Subset data to remove HH4
    SUBSET( TRANSFER_LABELS.out )

    // Recluster data
    CLUSTER_FULL( SUBSET.out )
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