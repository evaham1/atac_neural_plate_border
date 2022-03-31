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

include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { ARCHR_INTEGRATE } from "$baseDir/subworkflows/local/archr_integration"

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

    METADATA( params.sample_sheet )

    // add gtf to cellranger output so can add annotations
    METADATA.out // METADATA.out: [[meta], [cellranger_output]]
        .combine(ch_reference)
        .map{[it[0], it[1] + it[2]]}
        .set {ch_metadata} // ch_metadata: [[meta], [cellranger_output, gtf]]

    // ARCHR: run processing + clustering + filtering + gene scores on full data
    ARCHR_PROCESSING ( ch_metadata ) //output = archr_filtered_full

    // ARCHR: run clustering + gene scores on individual stages
    ARCHR_STAGE_PROCESSING ( ARCHR_PROCESSING.out.output )

    // ATAC: combine stage data and full data
    ARCHR_STAGE_PROCESSING.out.output
        .combine(ARCHR_PROCESSING.out.output)
        .view()
        .set {ch_atac}

    // RNA: read in data
    //.set ch_rna
    METADATA_RNA( params.rna_sample_sheet )
        .set {ch_rna}
        //.view()
   
    // combine ATAC and RNA data
    ch_atac
        .concat(ch_rna)
        .groupTuple(by:0)
        .view()
        .set {ch_integrate}

    // ARCHR: Integrate
    ARCHR_INTEGRATE(ch_integrate)
}

//[[sample_id:ss8], [[rds_files], [ss8_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:ss4], [[rds_files], [ss4_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:HH5], [[rds_files], [HH5_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:HH7], [[rds_files], [HH7_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:HH6], [[rds_files], [HH6_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:full], [[seurat_label_transfer.RDS]]]

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
