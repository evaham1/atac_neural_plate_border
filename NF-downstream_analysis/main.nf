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

// PREPROCESSING WORKFLOWS
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PREPROCESSING } from "$baseDir/subworkflows/local/Processing/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/Processing/Filtering"

// CLUSTERING WORKFLOW
include { CLUSTERING } from "$baseDir/subworkflows/local/archr_clustering_and_gex"

// INTEGRATION WORKFLOWS
include { METADATA as METADATA_ATAC } from "$baseDir/subworkflows/local/metadata"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { INTEGRATING } from "$baseDir/subworkflows/local/archr_integration"

// PEAK CALLING WORKFLOW
include { PEAK_CALLING } from "$baseDir/subworkflows/local/archr_peak_calling"

// PEAK EXPLORING WORKFLOWS
include { PEAK_EXPLORING } from "$baseDir/subworkflows/local/archr_peak_exploring"

// PARAMS
def skip_QC = params.skip_QC ? true : false

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

    if(!skip_QC){

        ///////////////////// PROCESSING //////////////////////////////
        METADATA( params.sample_sheet )    
        METADATA.out // METADATA.out: [[meta], [cellranger_output]]
            .combine(ch_reference)
            .map{[it[0], it[1] + it[2]]}
            .set {ch_metadata} // add gtf to cellranger output so can add annotations

        PREPROCESSING ( ch_metadata ) // create ArchR object

        FILTERING ( PREPROCESSING.out.output ) // iterative filtering
        ch_atac = FILTERING.out.output
        
    } else {
       
       METADATA_ATAC( params.atac_sample_sheet )
       ch_atac = METADATA_ATAC.out.metadata 
       // [[sample_id:HH5], [HH5_Save-ArchR]]
       //[[sample_id:HH6], [HH6_Save-ArchR]]
       //[[sample_id:HH7], [HH7_Save-ArchR]]
       //[[sample_id:ss4], [ss4_Save-ArchR]]
       //[[sample_id:ss8], [ss8_Save-ArchR]]
       //[[sample_id:FullData], [FullData_Save-ArchR]]

    }

    ///////////////////// CLUSTERING ////////////////////////////
    ///////////////////////////////////////////////////////////////
    CLUSTERING( ch_atac )

    ///////////////////// INTEGRATING //////////////////////////////
    ///////////////////////////////////////////////////////////////

    // RNA: read in data
    METADATA_RNA( params.rna_sample_sheet )
    //[[sample_id:HH5], [HH5_clustered_data.RDS]]
    //[[sample_id:HH6], [HH6_clustered_data.RDS]]
    //[[sample_id:HH7], [HH7_clustered_data.RDS]]
    //[[sample_id:ss4], [ss4_clustered_data.RDS]]
    //[[sample_id:ss8], [ss8_clustered_data.RDS]]
    //[[sample_id:FullData], [seurat_label_transfer.RDS]]
   
    // combine ATAC and RNA data
    CLUSTERING.out // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
        .concat( METADATA_RNA.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
        .groupTuple( by:0 ) //[ [sample_id:HH5], [ [rds_files], [HH5_splitstage_data/rds_files/HH5_clustered_data.RDS] ] ]
        .map{ [ it[0], [ it[1][0][3], it[1][1][0] ] ] }
        .view()
        .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

    // ARCHR: Integrate + filter out contaminating cells
    INTEGRATING( ch_integrate )
    // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]]]
    // [[ meta:full], ATAC]]
    // [ [[meta: HH5], ATAC] , [[meta: HH6], ATAC]]
    // [ [[meta: HH5], RNA], [[meta: HH6], RNA]]

    ///////////////////// PEAK CALLING ////////////////////////////
    ///////////////////////////////////////////////////////////////
    PEAK_CALLING( INTEGRATING.out.integrated_filtered )

    ///////////////////// EXPLORING //////////////////////////////
    ///////////////////////////////////////////////////////////////
    
    // IN PROCESS: peak exploring
    PEAK_EXPLORING( PEAK_CALLING.out )
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
