#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/NF-hichip-downstream
========================================================================================
    Github : 
    Website: 
    Slack  : 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

// 1) Create bins

// 2) Prep ValidPairs output from HiC pipeline and run with HiCDC+ to find significant interactions
include { METADATA as METADATA } from "$baseDir/subworkflows/local/metadata"
include {EDIT_VALIDPAIRS} from "$baseDir/modules/local/edit_ValidPairs/main"

// 2) Prep peak output and gtf and then intersect with bins to annotate them

// 3) Filter bins to pick out interactions of interest


// PARAMS


//
// SET CHANNELS
//

// set channel to NF-HiC output directory containing ValidPairs
// Channel
//     .value(params.)
//     .set{ch_reference}


//
// WORKFLOW: Run main nf-core/NF-hichip-downstream analysis pipeline
//

workflow A {

    // Read in ValidPairs data
    METADATA( params.sample_sheet_validpairs )

    // Edit ValidPairs data to add 'chr' to chromosome names
    EDIT_VALIDPAIRS ( METADATA.out.output )

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
    A ()
}


/*
========================================================================================
    THE END
========================================================================================
*/
