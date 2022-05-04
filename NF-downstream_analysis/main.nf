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
include { PREPROCESSING } from "$baseDir/subworkflows/local/1_processing/Preprocessing"
include { QC_STAGES as QC_NO_FITER } from "$baseDir/subworkflows/local/1_processing/Stage_processing"
include { QC_STAGES as QC_LOW } from "$baseDir/subworkflows/local/1_processing/Stage_processing"
include { QC_STAGES as QC_MED } from "$baseDir/subworkflows/local/1_processing/Stage_processing"
include { QC_STAGES as QC_HIGH } from "$baseDir/subworkflows/local/1_processing/Stage_processing"
include { FULL_PROCESSING as FULL_PROCESSING } from "$baseDir/subworkflows/local/1_processing/Full_processing"

// INTEGRATION WORKFLOWS
include { METADATA as METADATA_ATAC } from "$baseDir/subworkflows/local/metadata"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { INTEGRATING } from "$baseDir/subworkflows/local/archr_integration"


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
        METADATA.out.view()
        METADATA.out // METADATA.out: [[meta], [cellranger_output]]
            .combine(ch_reference)
            .map{[it[0], it[1] + it[2]]}
            .set {ch_metadata} // add gtf to cellranger output so can add annotations

        PREPROCESSING ( ch_metadata ) // create ArchR object

        //QC_NO_FILTER ( PREPROCESSING.out.output )
        //QC_LOW ( PREPROCESSING.out.output )
        QC_MED ( PREPROCESSING.out.output )
        //QC_HIGH ( PREPROCESSING.out.output )

        ch_combined = QC_MED.out.output // Collect rds files from all stages
            .concat(PREPROCESSING.out.output)
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            .collect()
            .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]
        FULL_PROCESSING ( ch_combined ) // filter full data
        ///////////////////////////////////////////////////////////////

        ch_atac = QC_MED.out.output // Collect rds files from all stages
            .concat(FULL_PROCESSING.out.output)
            .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //[ [[meta: HH5], ATAC.rds] , [[meta: HH6], ATAC.rds], [[meta: FullData], ATAC.rds]]
            .view()
        
    } else {
       
       METADATA_ATAC( params.atac_sample_sheet )
       ch_atac = METADATA_ATAC.out.metadata //[ [[meta: HH5], ATAC.rds] , [[meta: HH6], ATAC.rds], [[meta: FullData], ATAC.rds]]

    }

    ///////////////////// INTEGRATING //////////////////////////////
    ///////////////////////////////////////////////////////////////

    ch_atac.view() 
    // RNA: read in data
    METADATA_RNA( params.rna_sample_sheet )
   
    // combine ATAC and RNA data
    ch_atac
        .concat( METADATA_RNA.out.metadata )
        .groupTuple( by:0 )
        // .map{[it[0], it[[1]] + it[2]]}
        .map{ [ it[0], [it[1][0], it[1][1][0]] ] }
        .view()
        .set {ch_integrate}

    // ARCHR: Integrate
    INTEGRATING( ch_integrate )
    // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]]]
    // [[ meta:full], ATAC]]
    // [ [[meta: HH5], ATAC] , [[meta: HH6], ATAC]]
    // [ [[meta: HH5], RNA], [[meta: HH6], RNA]]

    // INTEGRATING.out.archr_integrated_full.view()
    
    // ///////////////////// PEAK CALLING ////////////////////////////
    // ///////////////////////////////////////////////////////////////
    
    // PEAK_CALLING( INTEGRATING.out.archr_integrated_full )
    
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
