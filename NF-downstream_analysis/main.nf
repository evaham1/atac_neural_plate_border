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
include { PREPROCESSING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Filtering"

// PROCESSING WORKFLOWS
include { CLUSTERING as CLUSTERING_WITH_CONTAM } from "$baseDir/subworkflows/local/PROCESSING/archr_clustering_and_gex"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { INTEGRATING } from "$baseDir/subworkflows/local/PROCESSING/archr_integration"
include { PEAK_CALLING } from "$baseDir/subworkflows/local/PROCESSING/archr_peak_calling"

include { TRANSFER_LABELS } from "$baseDir/subworkflows/local/PROCESSING/archr_transfer_labels"

// DOWNSTREAM PROCESSING WORKFLOWS
include { COMPARE_VARIABILITY } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_compare_variability"
include { NPB_SUBSET } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_npb_subset"

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
    // creates transfer labels object

    if(!skip_processing){

        //ch_upstream_processed.view()

        // Cluster QC'd atac cells
        CLUSTERING_WITH_CONTAM( ch_upstream_processed )
        //CLUSTERING_WITH_CONTAM.out.output.view()

        // Extract the stages to run integration on them
        CLUSTERING_WITH_CONTAM.out.output
            .filter{ meta, data -> meta.sample_id != 'FullData'} // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()]}
            .set{ ch_stages } // [ [sample_id:HH5], [HH5-ArchR] ]
             

        // read in RNA data
        METADATA_RNA( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                             // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                             // etc
   

        // combine ATAC and RNA data
        ch_stages
            .concat( METADATA_RNA.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
            .groupTuple( by:0 )
            .map{ meta, data -> [meta, [data[0][0], data[1][0]]]}
            //.view()
            .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

        

        // ARCHR: Integrate + filter out contaminating cells
        INTEGRATING( ch_integrate )  // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]], etc]

        // Call peaks on resulting data (stages + full filtered for contamination)
        //PEAK_CALLING( INTEGRATING.out.integrated_filtered ) //TEMP COMMENTED OUT

        /////////////// Transfer labels from stages onto full data  //////////////////////////

        // extract the full data
        CLUSTERING_WITH_CONTAM.out.view()
        CLUSTERING_WITH_CONTAM.out
            .filter{ meta, data -> meta.sample_id == 'FullData'}
            .set{ ch_fulldata_clustered }
            //.view()

        INTEGRATING.out.integrated_filtered
            .concat( ch_fulldata_clustered )
            .map{ meta, data -> [data.findAll{it =~ /rds_files/}[0].listFiles()] } //removes all metadata and list files in rds_files
            .collect()
            .map{data -> [[sample_id:'transfer_labels'], [data]] }
            .view()
            .set{ ch_transfer_labels_input }



        // and combine to do transfer labels

        // ch_fulldata_clustered
        //     .concat{ ch_stages_integrated }
        //     .map{ meta, data -> [data.findAll{it =~ /rds_files/}[0].listFiles()] } //removes all metadata and list files in rds_files
        //     .collect() //channel of length 6 turns into channel length 1
        //     .map{ data -> [[sample_id:'transfer_labels'], [data]] }
        //     .set{ ch_transfer_labels_input }
        // //TRANSFER_LABELS( ch_transfer_labels_input )

        //ch_processed = INTEGRATING.out.integrated_filtered //TEMP
        //ch_processed = PEAK_CALLING.out
        //ch_processed_transfer_labels = TRANSFER_LABELS.out

    } //else {
       
    //    METADATA_PROCESSED( params.processed_sample_sheet )
    //    // ADD TRANSFER LABELS OBJECT TO THIS SAMPLE SHEET
    //    ch_processed = METADATA_PROCESSED.out.metadata                       // [[sample_id:HH5], [HH5_Save-ArchR]]
    //                                                                         // [[sample_id:HH6], [HH6_Save-ArchR]]
    //                                                                         // etc

    // }


    ///////////////////////////////////////////////////////////////
    ///////////////////// DOWNSTREAM PROCESSING ///////////////////
    ///////////////////////////////////////////////////////////////
    // comparing stages
    // making transfer_labels full object and working on that
    // WORK IN PROGRESS
    
    // IN PROGRESS: compare variability of clusters between stages
    // currently just uses differential peak tests, would be better to measure in another way
    //COMPARE_VARIABILITY( ch_processed )

    // IN PROGRESS: combine individual stages integrated objects into full transferlabels object and look for enhancers/study peaks over time
    //TRANSFER_LABELS( ch_processed )
    //Unexpected error [StackOverflowError]
    //-- Check script 'subworkflows/local/DOWNSTREAM_PROCESSING/archr_transfer_labels.nf' at line: 48 or see '.nextflow.log' file for more details

    // IN PROGRESS: subset out NPB subset from transfer labels object and focus on that
    //NPB_SUBSET( TRANSFER_LABELS? )

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
