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

include { METADATA_WF } from "$baseDir/subworkflows/local/metadata"

include { PROCESSING_ARCHR_WF } from "$baseDir/subworkflows/local/archr_processing"
include { STAGE_PROCESSING_ARCHR_WF } from "$baseDir/subworkflows/local/archr_stage_processing"

include { METADATA_WF as METADATA_WF_RNA } from "$baseDir/subworkflows/local/metadata"
include { INTEGRATION_ARCHR_WF } from "$baseDir/subworkflows/local/archr_integration"

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
workflow A {

    ///////////////////// PROCESSING //////////////////////////////
    ///////////////////////////////////////////////////////////////
    
    METADATA_WF( params.sample_sheet )

    // add gtf to cellranger output so can add annotations
    METADATA_WF.out // METADATA.out: [[meta], [cellranger_output]]
        .combine(ch_reference)
        .map{[it[0], it[1] + it[2]]}
        .set {ch_metadata} // ch_metadata: [[meta], [cellranger_output, gtf]]

    // ARCHR: run processing + clustering + filtering + gene scores on full data
    PROCESSING_ARCHR_WF ( ch_metadata ) //output = archr_filtered_full

    // ARCHR: run clustering + gene scores on individual stages
    STAGE_PROCESSING_ARCHR_WF ( PROCESSING_ARCHR_WF.out.output )

    // ATAC: add together stage data and full data
    STAGE_PROCESSING_ARCHR_WF.out.output
        .concat( PROCESSING_ARCHR_WF.out.output )
        //.view()
        .set {ch_atac}

    ///////////////////// INTEGRATING //////////////////////////////
    ///////////////////////////////////////////////////////////////

    // RNA: read in data
    METADATA_WF_RNA( params.rna_sample_sheet )
   
    // combine ATAC and RNA data
    ch_atac
        .concat( METADATA_WF_RNA.out.metadata )
        .groupTuple( by:0 )
        // .map{[it[0], it[[1]] + it[2]]}
        .map{ [ it[0], [it[1][0], it[1][1][0]] ] }
        //.view()
        .set {ch_integrate}

    // ARCHR: Integrate
    INTEGRATION_ARCHR_WF( ch_integrate )

    INTEGRATION_ARCHR_WF.out.archr_integrated_full.view()
    
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
