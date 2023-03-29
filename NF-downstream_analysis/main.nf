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
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_split_stages.R", checkIfExists: true) )

//CALCULATING SEACELL WFs
include { SEACELLS_ATAC_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_ATAC_WF"

include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { SEACELLS_RNA_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_RNA_WF"

//INTEGRATING SEACELLS
include {PYTHON as INTEGRATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/Integration/SEACells_integration.py", checkIfExists: true) )


// DOWNSTREAM PROCESSING WORKFLOWS

// include { CLUSTER_PEAKS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/cluster_peaks_WF"
// include { FIND_ENHANCERS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/find_enhancers_WF"
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

// Set channel for binary knowledge matrix for cell state classification
Channel
    .value("$baseDir/binary_knowledge_matrix_contam.csv")
    .set{ch_binary_knowledge_matrix}


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

        ch_metadata.view()

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
    // calls peaks on full data object to obtain consensus peak set and peak counts
    // calculates metacells for RNA and ATAC stages
    // run integration between RNA and ATAC metacells

    if(!skip_processing){

        ///////     Call peaks on full data      ///////

        // Extract just the full data object
        ch_upstream_processed
            .filter{ meta, data -> meta.sample_id == 'FullData'}
            .set{ ch_full }

        // Cluster and call peaks on full data
        //ch_full.view() //[[sample_id:FullData], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/output/NF-downstream_analysis/Upstream_processing/FILTERING/FullData/rds_files/FullData_Save-ArchR]]
        CLUSTER( ch_full )
        PEAK_CALL( CLUSTER.out )

        ///////     Run Metacells on stage data      ///////

        // Split full data into stages
        SPLIT_STAGES( PEAK_CALL.out )
        SPLIT_STAGES.out //[[meta], [plots, rds_files]]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .flatMap {it[1][0].listFiles()}
            .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
            .set { ch_split_stage }

        // Run Metacells on ATAC stages
        ch_split_stage.view()
        SEACELLS_ATAC_WF( ch_split_stage, ch_binary_knowledge_matrix )
             
        //read in RNA data (stages only)
        METADATA_RNA( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                                // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                                // etc
        // Run Metacells on RNA stages
        SEACELLS_RNA_WF( METADATA_RNA.out.metadata, ch_binary_knowledge_matrix )

        ///////     Integrate SEACells      ///////

        //SEACELLS_RNA_WF.out.seacells_anndata_processed_classified.view()
        // [[sample_id:HH6], [c8/0355d43ba613bc2159e10758c26572/plots, c8/0355d43ba613bc2159e10758c26572/rds_files]]
        // [[sample_id:HH5], [61/920e02828bc1d30d01745280da08cb/plots, 61/920e02828bc1d30d01745280da08cb/rds_files]]
        // [[sample_id:ss8], [64/bc0c8d3f045086356bc00a7ebde6ff/plots, 64/bc0c8d3f045086356bc00a7ebde6ff/rds_files]]
        // [[sample_id:ss4], [fa/1d64469d41d5a69f74d4d706bf5681/plots, fa/1d64469d41d5a69f74d4d706bf5681/rds_files]]
        // [[sample_id:HH7], [d0/6b148d2543b6b195b8ee2088fad0b7/plots, d0/6b148d2543b6b195b8ee2088fad0b7/rds_files]]

        SEACELLS_RNA_WF.out.seacells_anndata_processed_classified
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()[0]] }
            .set {ch_RNA}

        // ch_RNA.view()
        // [[sample_id:HH6], c8/0355d43ba613bc2159e10758c26572/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:HH5], 61/920e02828bc1d30d01745280da08cb/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:HH7], d0/6b148d2543b6b195b8ee2088fad0b7/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:ss4], fa/1d64469d41d5a69f74d4d706bf5681/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:ss8], 64/bc0c8d3f045086356bc00a7ebde6ff/rds_files/AnnData_RNA.h5ad]
        
        SEACELLS_ATAC_WF.out.seacells_anndata_processed_classified
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()[0]] }
            .set {ch_ATAC}
        // ch_ATAC is identical to ch_RNA except objects are called AnnData_ATAC.h5ad


        ch_RNA
            .concat( ch_ATAC )
            .groupTuple( by:0 )
            .view() //[[sample_id:ss4], [fa/1d64469d41d5a69f74d4d706bf5681/rds_files/AnnData_RNA.h5ad, 4a/72711280da408db03e302a5926456a/rds_files/AnnData_ATAC.h5ad]]
            .set {ch_seacells_to_integrate}
        INTEGRATE_SEACELLS( ch_seacells_to_integrate )






        // /////////////// Transfer labels from integrated stages onto non-integrated full data  //////////////////////////

        // // extract the full data
        // CLUSTER.out
        //     .filter{ meta, data -> meta.sample_id == 'FullData'}
        //     //.view() //[[sample_id:FullData], [./rds_files]]
        //     .set{ ch_fulldata_clustered }

        // // combine clustered full data with integrated stage data into one channel
        // INTEGRATING.out.integrated_filtered
        //     .concat( ch_fulldata_clustered )
        //     .map{ meta, data -> [data.findAll{it =~ /rds_files/}[0].listFiles()[0]] } //removes all metadata and list files in rds_files
        //     //.view()
        //             //[HH5_Save-ArchR]
        //             //[HH7_Save-ArchR]
        //             //[ss4_Save-ArchR]
        //             //[ss8_Save-ArchR]
        //             //[HH6_Save-ArchR]
        //             //[FullData_Save-ArchR]
        //     .collect()
        //     .map{data -> [[sample_id:'TransferLabels'], data] }
        //     //.view() //[[sample_id:TransferLabels], [[HH5_Save-ArchR, HH7_Save-ArchR, ss4_Save-ArchR, ss8_Save-ArchR, HH6_Save-ArchR, FullData_Save-ArchR]]]
        //     .set{ ch_transfer_labels_input }

        // // transfers labels to full object, clusters and call peaks on stage_clusters of transferlabels object
        // TRANSFER_LABELS( ch_transfer_labels_input )

        // /////////////// Create output channel  //////////////////////////
        // CLUSTER.out
        //     .concat( TRANSFER_LABELS.out.transfer_label_peaks )
        //     //.view()
        //     .set{ ch_processed }

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
    // ch_processed
    //         .filter{ meta, data -> meta.sample_id == 'TransferLabels'}
    //         //.view() //[[sample_id:TransferLabels], [./TransferLabels_Save-ArchR]]
    //         .set{ ch_TL }

    // // Subworkflow to create metacells
    // SEACELLS_WF( ch_TL )

    // // Subworkflow to cluster peaks using metacells
    // CLUSTER_PEAKS_WF( CALCULATE_SEACELLS.out.seacells_output_combined )

    // // Subworkflow to identify NC and PPR specific enhancers ~ maybe don't need if can find these in the peak modules?
    // FIND_ENHANCERS_WF( ch_TL )

    
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
