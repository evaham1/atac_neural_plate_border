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

// METADATA WORKFLOWS FOR CHANNEL SWITCHES
include { METADATA as METADATA_UPSTREAM_PROCESSED } from "$baseDir/subworkflows/local/metadata"
include { METADATA as METADATA_PROCESSED } from "$baseDir/subworkflows/local/metadata"

// UPSTREAM PROCESSING WORKFLOWS
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PREPROCESSING } from "$baseDir/subworkflows/local/Processing/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/Processing/Filtering"

// PROCESSING WORKFLOWS
include { CLUSTERING as CLUSTERING_WITH_CONTAM } from "$baseDir/subworkflows/local/archr_clustering_and_gex"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { INTEGRATING } from "$baseDir/subworkflows/local/archr_integration"
include { PEAK_CALLING } from "$baseDir/subworkflows/local/archr_peak_calling"

// DOWNSTREAM PROCESSING WORKFLOWS
include { PEAK_EXPLORING } from "$baseDir/subworkflows/local/archr_peak_exploring"

// PARAMS
def skip_upstream_processing = params.skip_upstream_processing ? true : false
def skip_processing = params.skip_processing ? true : false

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

    ///////////////////////////////////////////////////////////////
    ///////////////////// UPSTREAM PROCESSING /////////////////////
    ///////////////////////////////////////////////////////////////
    // sets up ArchR object
    // QC and filtering

workflow A {

    if(!skip_upstream_processing){

        METADATA( params.sample_sheet )    
        METADATA.out // METADATA.out: [[meta], [cellranger_output]]
            .combine(ch_reference)
            .map{[it[0], it[1] + it[2]]}
            .set {ch_metadata} // add gtf to cellranger output so can add annotations

        PREPROCESSING ( ch_metadata ) // create ArchR object

        FILTERING ( PREPROCESSING.out.output ) // iterative filtering
        ch_upstream_processed = FILTERING.out.output
        
    } else {
       
       METADATA_UPSTREAM_PROCESSED( params.upstream_processed_sample_sheet )
       ch_upstream_processed = METADATA_UPSTREAM_PROCESSED.out.metadata     // [[sample_id:HH5], [HH5_Save-ArchR]]
                                                                            // [[sample_id:HH6], [HH6_Save-ArchR]]
                                                                            // etc

    }

    ///////////////////////////////////////////////////////////////
    /////////////////////    PROCESSING      //////////////////////
    ///////////////////////////////////////////////////////////////
    // integrates with scRNA, filters out contam
    // clusters
    // calls peaks

    if(!skip_processing){

        // Cluster QC'd atac cells
        CLUSTERING_WITH_CONTAM( ch_upstream_processed )

        // read in RNA data
        METADATA_RNA( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                            // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                            // etc
   
        // combine ATAC and RNA data
        CLUSTERING_WITH_CONTAM.out // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .concat( METADATA_RNA.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
            .groupTuple( by:0 ) //[ [sample_id:HH5], [ [rds_files], [HH5_splitstage_data/rds_files/HH5_clustered_data.RDS] ] ]
            .map{ [ it[0], [ it[1][0][3], it[1][1][0] ] ] }
            .view()
            .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

        // ARCHR: Integrate + filter out contaminating cells
        INTEGRATING( ch_integrate )  // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]], etc]

        // Call peaks on resulting data
        PEAK_CALLING( INTEGRATING.out.integrated_filtered )

        ch_processed = PEAK_CALLING.out

    } else {
       
       METADATA_PROCESSED( params.processed_sample_sheet )
       ch_processed = METADATA_PROCESSED.out.metadata                       // [[sample_id:HH5], [HH5_Save-ArchR]]
                                                                            // [[sample_id:HH6], [HH6_Save-ArchR]]
                                                                            // etc

    }

    ch_processed.view()


    ///////////////////////////////////////////////////////////////
    ///////////////////// DOWNSTREAM PROCESSING ///////////////////
    ///////////////////////////////////////////////////////////////
    // comparing stages
    // making transfer_labels full object and working on that
    // WORK IN PROGRESS
    
    // IN PROCESS: peak exploring
    PEAK_EXPLORING( ch_processed )
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
