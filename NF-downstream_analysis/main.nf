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

include { ARCHR_PROCESSING } from "$baseDir/subworkflows/local/archr_processing"
include { ARCHR_STAGE_PROCESSING } from "$baseDir/subworkflows/local/archr_stage_processing"

include { ARCHR_INTEGRATE } from "$baseDir/subworkflows/local/archr_integration"

//
// SET CHANNELS
//

// set channel to reference folder containing fasta and gtf
Channel
    .value(params.reference)
    .set{ch_reference}

// set channel to rna RDS object - will need to adjust this so has contaminating clusters
Channel
    .value(params.seurat_RNA)
    .set{ch_rna} //ch_rna: seurat.RDS


//
// WORKFLOW: Run main nf-core/downstream analysis pipeline
//
workflow NFCORE_DOWNSTREAM {

    METADATA( params.sample_sheet )

    // add gtf to cellranger output so can add annotations
    METADATA.out // METADATA.out: [[meta], [cellranger_output]]
        .combine(ch_reference)
        .map{[it[0], it[1] + it[2]]}
        .set {ch_metadata} // ch_metadata: [[meta], [cellranger_output, gtf]]

    // ARCHR: run processing + clustering + filtering + gene scores on full data
    ARCHR_PROCESSING ( ch_metadata ) //output = archr_filtered_full

    // ARCHR: run clustering + gene scores on individual stages
    ARCHR_STAGE_PROCESSING ( ARCHR_PROCESSING.out.archr_filtered_full ) //output = archr_filtered_stages

    // add full ATAC data to list of ATAC stage data
    // ARCHR_PROCESSING.out.archr_filtered_full
    //     .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]}
    //     .view()
    //     .set { ch_atac_full }

    // ARCHR_STAGE_PROCESSING.out.atac_stage_merged
    //     .combine(ch_atac_full) // Combine with full dataset
    //     .map{ row -> [row[0], [row[1], row[2]]]}
    //     .view()
    //     .set { ch_atac_integrate } 

    // // read into integration 
    // ARCHR_INTEGRATE( ch_atac_integrate )
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
