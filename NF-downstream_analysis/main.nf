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
include { INTEGRATE_SPLIT_PROCESS } from "$baseDir/subworkflows/local/integrate_split_process"
include {R as INTEGRATE_RNA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/4_integrate_rna.R", checkIfExists: true) )

// set channel to reference gtf
Channel
    .value(params.gtf)
    .set{ch_gtf}

// set channel to rna RDS object
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
        .combine(ch_gtf)
        .map{[it[0], it[1] + it[2]]}
        .set {ch_metadata} // ch_metadata: [[meta], [cellranger_output, gtf]]
    
    // run preprocessing, filtering and predicted gex
    PROCESSING ( ch_metadata )

    // strip metadata from outputs
    PROCESSING.out.signac_predicted_gex
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .set {ch_atac} // ch_atac: seurat_GeneActivity.RDS
    METADATA.out
        .map{it[1]}
        .set {ch_cellranger} // ch_cellranger: [cellranger_output]
    
    // run rna integration on individual stages
    INTEGRATE_SPLIT_PROCESS( ch_atac , ch_rna, ch_cellranger )

    INTEGRATE_SPLIT_PROCESS.out.signac_integrated.view()
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
