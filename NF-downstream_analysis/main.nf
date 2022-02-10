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
    GENOME PARAMETER VALUES
========================================================================================
*/

//params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

//WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PROCESSING } from "$baseDir/subworkflows/local/processing"

// set channel to reference gtf
Channel
    .value(params.gtf)
    .set{ch_gtf}

// set channel to rna RDS object
Channel
    .value(params.transfer_labels)
    .set{ch_transfer_labels}

//
// WORKFLOW: Run main nf-core/downstream analysis pipeline
//
workflow NFCORE_DOWNSTREAM {
    //ch_gtf.view()

    METADATA( params.sample_sheet )

    //METADATA.out.view()

    METADATA.out
        .combine(ch_gtf)
        .map{[it[0], it[1] + it[2]]}
        .view()
        .set {ch_metadata}
    
    PROCESSING (ch_metadata )
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
