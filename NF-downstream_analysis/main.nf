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

// set inputs 
// input = [
//         [ id:'test'],
//         "$baseDir/../output/NF-luslab_sc_multiomic/hh7_1_cellranger_atac/outs/*"
//     ]

//
// WORKFLOW: Run main nf-core/downstream analysis pipeline
//
workflow NFCORE_DOWNSTREAM {
    METADATA( 'samplesheet.csv' )

    // filter metadata outputs to only keep filtered_peak_bc_matrix.h5, singlecell.csv, fragments.tsv.gz
    // METADATA.out
    //     .filter{ it[0].sample_id == 'NF-scATACseq_alignment_out' }
    //     .map {[it[0], it[1].collect{ file(it+"/outs/filtered_peak_bc_matrix", checkIfExists: true) }]}
    //     .set {ch_scATACseq_peaks}
    
    PROCESSING ( METADATA.out )
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
