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

include { WF_METADATA } from "$baseDir/subworkflows/local/metadata"

include { WF_ARCHR_PROCESSING } from "$baseDir/subworkflows/local/archr_processing"
include { WF_ARCHR_STAGE_PROCESSING } from "$baseDir/subworkflows/local/archr_stage_processing"

include { WF_METADATA as WF_METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { WF_ARCHR_INTEGRATION } from "$baseDir/subworkflows/local/archr_integration"

//
// SET CHANNELS
//

// set channel to reference folder containing fasta and gtf
Channel
    .value(params.reference)
    .set{ch_reference}


//
// WORKFLOW: Run main nf-core/downstream analysis pipeline
//
workflow NFCORE_DOWNSTREAM {

    ///////////////////// PROCESSING //////////////////////////////
    ///////////////////////////////////////////////////////////////
    
    WF_METADATA( params.sample_sheet )

    // add gtf to cellranger output so can add annotations
    WF_METADATA.out // METADATA.out: [[meta], [cellranger_output]]
        .combine(ch_reference)
        .map{[it[0], it[1] + it[2]]}
        .set {ch_metadata} // ch_metadata: [[meta], [cellranger_output, gtf]]

    // ARCHR: run processing + clustering + filtering + gene scores on full data
    WF_ARCHR_PROCESSING ( ch_metadata ) //output = archr_filtered_full

    // ARCHR: run clustering + gene scores on individual stages
    WF_ARCHR_STAGE_PROCESSING ( WF_ARCHR_PROCESSING.out.output )

    // ATAC: add together stage data and full data
    WF_ARCHR_STAGE_PROCESSING.out.output
        .concat( WF_ARCHR_PROCESSING.out.output )
        //.view()
        .set {ch_atac}

    ///////////////////// INTEGRATING //////////////////////////////
    ///////////////////////////////////////////////////////////////

    // RNA: read in data
    WF_METADATA_RNA( params.rna_sample_sheet )
   
    // combine ATAC and RNA data
    ch_atac
        .concat( WF_METADATA_RNA.out.metadata )
        .groupTuple( by:0 )
        // .map{[it[0], it[[1]] + it[2]]}
        .map{ [ it[0], [it[1][0], it[1][1][0]] ] }
        //.view()
        .set {ch_integrate}

    // ARCHR: Integrate
    WF_ARCHR_INTEGRATION( ch_integrate )
    
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
