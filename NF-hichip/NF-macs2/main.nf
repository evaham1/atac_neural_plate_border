#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/macs2
========================================================================================
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

// LOAD MODULES
include { METADATA } from "./subworkflows/local/metadata"
include { MACS2_CALLPEAK           } from './modules/nf-core/macs2/callpeak/main'

// SET CHANNELS
Channel
    .value(params.macs_gsize)
    .set{macs_gsize}


//ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
//macs_gsize                        // integer: value for --macs_gsize parameter

//
// WORKFLOW: Run macs2
//

workflow NF_PEAK_CALL {

    METADATA( params.sample_sheet )

    MACS2_CALLPEAK (METADATA.out, macs_gsize)

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
    NF_PEAK_CALL ()
}


/*
========================================================================================
    THE END
========================================================================================
*/
