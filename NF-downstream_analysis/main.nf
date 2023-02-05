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

// PROCESSING WORKFLOWS AND MODULES
//MODULES
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )

//WORKFLOWS
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { INTEGRATING } from "$baseDir/subworkflows/local/PROCESSING/archr_integration"
include { TRANSFER_LABELS } from "$baseDir/subworkflows/local/PROCESSING/archr_transfer_labels"


// DOWNSTREAM PROCESSING WORKFLOWS
include { CALCULATE_SEACELLS } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/calculate_seacells"
//include { COMPARE_VARIABILITY } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_compare_variability"
//include { NPB_SUBSET } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_npb_subset"

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

        ///
        ch_upstream_processed = FILTERING.out.output
        //ch_upstream_processed.view()
            //[[sample_id:ss4], [ss4_Save-ArchR]]
            //[[sample_id:HH7], [HH7_Save-ArchR]]
            //[[sample_id:ss8], [ss8_Save-ArchR]]
            //[[sample_id:HH5], [HH5_Save-ArchR]]
            //[[sample_id:HH6], [HH6_Save-ArchR]]
            //[[sample_id:FullData], [FullData_Save-ArchR]]
        ///
        
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

        /////////////// Cluster individual stages and full data  //////////////////////////

        // Cluster stages + full data
        CLUSTER( ch_upstream_processed )

        /////////////// Integrate stages data with RNA stages data  //////////////////////////

        // Extract the stages to run integration on them
        CLUSTER.out
            .filter{ meta, data -> meta.sample_id != 'FullData'} // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()]}
            .set{ ch_stages } // [ [sample_id:HH5], [HH5-ArchR] ]
             
        // read in RNA data (stages only)
        METADATA_RNA( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                                // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                                // etc

        // combine ATAC and RNA data (stages only)
        ch_stages
            .concat( METADATA_RNA.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
            .groupTuple( by:0 )
            .map{ meta, data -> [meta, [data[0][0], data[1][0]]]}
            //.view()
            .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

        // Integrate + filter out contaminating cells (stages only)
        INTEGRATING( ch_integrate )  // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]], etc]

        /////////////// Call peaks on integrated, contam filtered stages data  //////////////////////////

        // Call peaks on resulting data (stages only)
        PEAK_CALL( INTEGRATING.out.integrated_filtered )

        /////////////// Transfer labels from integrated stages onto non-integrated full data  //////////////////////////

        // extract the full data
        CLUSTER.out
            .filter{ meta, data -> meta.sample_id == 'FullData'}
            //.view() //[[sample_id:FullData], [./rds_files]]
            .set{ ch_fulldata_clustered }

        // combine clustered full data with integrated stage data into one channel
        INTEGRATING.out.integrated_filtered
            .concat( ch_fulldata_clustered )
            .map{ meta, data -> [data.findAll{it =~ /rds_files/}[0].listFiles()[0]] } //removes all metadata and list files in rds_files
            //.view()
                    //[HH5_Save-ArchR]
                    //[HH7_Save-ArchR]
                    //[ss4_Save-ArchR]
                    //[ss8_Save-ArchR]
                    //[HH6_Save-ArchR]
                    //[FullData_Save-ArchR]
            .collect()
            .map{data -> [[sample_id:'TransferLabels'], data] }
            //.view() //[[sample_id:TransferLabels], [[HH5_Save-ArchR, HH7_Save-ArchR, ss4_Save-ArchR, ss8_Save-ArchR, HH6_Save-ArchR, FullData_Save-ArchR]]]
            .set{ ch_transfer_labels_input }

        // transfers labels to full object, clusters and call peaks on stage_clusters of transferlabels object
        TRANSFER_LABELS( ch_transfer_labels_input )

        /////////////// Create output channel  //////////////////////////
        CLUSTER.out
            .concat( TRANSFER_LABELS.out.transfer_label_peaks )
            //.view()
            .set{ ch_processed }

    } else {
       
       METADATA_PROCESSED( params.processed_sample_sheet )
       // !! NEED TO ADD TRANSFER LABELS OBJECT TO THIS SAMPLE SHEET
       ch_processed = METADATA_PROCESSED.out.metadata                       
                                                                            //[[sample_id:HH5], [/camp/home/hamrude/scratch/atac_neural_plate_border/output/NF-downstream_analysis/Processing/HH5/7_peak_call/rds_files/HH5_Save-ArchR]]
                                                                            //[[sample_id:HH6], [/camp/home/hamrude/scratch/atac_neural_plate_border/output/NF-downstream_analysis/Processing/HH6/7_peak_call/rds_files/HH6_Save-ArchR]]
                                                                            //[[sample_id:HH7], [/camp/home/hamrude/scratch/atac_neural_plate_border/output/NF-downstream_analysis/Processing/HH7/7_peak_call/rds_files/HH7_Save-ArchR]]
                                                                            //[[sample_id:ss4], [/camp/home/hamrude/scratch/atac_neural_plate_border/output/NF-downstream_analysis/Processing/ss4/7_peak_call/rds_files/ss4_Save-ArchR]]
                                                                            //[[sample_id:ss8], [/camp/home/hamrude/scratch/atac_neural_plate_border/output/NF-downstream_analysis/Processing/ss8/7_peak_call/rds_files/ss8_Save-ArchR]]
                                                                            //[[sample_id:TransferLabels], [/camp/home/hamrude/scratch/atac_neural_plate_border/output/NF-downstream_analysis/Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]

    }


    ///////////////////////////////////////////////////////////////
    ///////////////////// DOWNSTREAM PROCESSING ///////////////////
    ///////////////////////////////////////////////////////////////


    //Extract just TransferLabels object from ch_processed
    ch_processed
            .filter{ meta, data -> meta.sample_id == 'TransferLabels'}
            //.view() //[[sample_id:TransferLabels], [./TransferLabels_Save-ArchR]]
            .set{ ch_TL }

    // Subworkflow to create metacells
    CALCULATE_SEACELLS( ch_TL )

    // Subworkflow to cluster peaks using metacells

    
    // IN PROGRESS: compare variability of clusters between stages
    // currently just uses differential peak tests, would be better to measure in another way
    //COMPARE_VARIABILITY( ch_processed )

    // IN PROGRESS: check validity of cell state labels
    // can try reclustering on differet peak sets and see if get more distinct label separation
    // script in integration/compare_clusters_and_labels
    // predict enhancers using diff accessibility and validate them
    // findenhancers script

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
