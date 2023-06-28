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
include { METADATA as METADATA_PEAKCALL_PROCESSED } from "$baseDir/subworkflows/local/metadata"
include { METADATA as METADATA_SINGLECELL_PROCESSED } from "$baseDir/subworkflows/local/metadata"

// UPSTREAM PROCESSING
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PREPROCESSING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Filtering"

// PROCESSING
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_peak_calling.R", checkIfExists: true) )
include {R as SPLIT_STAGES_PROCESSED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_split_stages.R", checkIfExists: true) )

// SINGLE CELL PROCESSING
include { METADATA as METADATA_RNA_SC } from "$baseDir/subworkflows/local/metadata"
include { ARCHR_STAGE_PEAKS_WF } from "$baseDir/subworkflows/local/PROCESSING/archr_stage_peaks_WF"
include { ARCHR_INTEGRATING_WF } from "$baseDir/subworkflows/local/PROCESSING/archr_integration_WF"
include { ARCHR_STAGE_DIFF_PEAKS_WF } from "$baseDir/subworkflows/local/PROCESSING/archr_stage_diff_peaks_WF"

// METACELL PROCESSING
include { SEACELLS_ATAC_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_ATAC_WF"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { SEACELLS_RNA_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_RNA_WF"
include { SEACELLS_INTEGRATING_WF } from "$baseDir/subworkflows/local/PROCESSING/SEACells_integration_WF"

//PEAK CLUSTERING
include { CLUSTER_PEAKS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/cluster_peaks_WF"

// MISC
include {R as MAKE_TXDB} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/gtf_to_txdb.R", checkIfExists: true) )
include {EXTRACT_EXONS} from "$baseDir/modules/local/extract_exons/main"

// DOWNSTREAM PROCESSING WORKFLOWS

// 
// include { FIND_ENHANCERS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/find_enhancers_WF"
//include { COMPARE_VARIABILITY } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_compare_variability"
//include { NPB_SUBSET } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_npb_subset"

// PARAMS FOR PIPELINE SWITCHING
def skip_upstream_processing = params.skip_upstream_processing ? true : false
def skip_peakcall_processing = params.skip_peakcall_processing ? true : false
def skip_metacell_processing = params.skip_metacell_processing ? true : false
def skip_singlecell_processing = params.skip_singlecell_processing ? true : false


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

// set channel to just gtf file
Channel
    .value(params.gtf)
    .set{ch_gtf}


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

        METADATA( params.aligned_sample_sheet )    
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


    //////////////////////////////////////////////////////////////////////////////////
    /////////////////////      GENERATE CONSENSUS PEAK SET      //////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    // calls peaks on full data object to obtain consensus peak set, then splits full data into stages

    if(!skip_peakcall_processing){

        // Extract just the full data object
        ch_upstream_processed
            .filter{ meta, data -> meta.sample_id == 'FullData'}
            .set{ ch_full }

        // Cluster and call peaks on full data
        //ch_full.view() //[[sample_id:FullData], [Upstream_processing/FILTERING/FullData/rds_files/FullData_Save-ArchR]]
        CLUSTER( ch_full )
        PEAK_CALL( CLUSTER.out )

        // Split full data into stages
        SPLIT_STAGES_PROCESSED( PEAK_CALL.out )
        SPLIT_STAGES_PROCESSED.out //[[meta], [plots, rds_files]]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .flatMap {it[1][0].listFiles()}
            .map { row -> [[sample_id:dummy], row] }
            .set { ch_peakcall_processed }

        SPLIT_STAGES_PROCESSED.out.view()

        // ch_peakcall_processed.view()
        // [[sample_id:HH6], rds_files/HH6_Save-ArchR]
        // [[sample_id:HH7], rds_files/HH7_Save-ArchR]
        // [[sample_id:ss4], rds_files/ss4_Save-ArchR]
        // [[sample_id:ss8], rds_files/ss8_Save-ArchR]
        // [[sample_id:HH5], rds_files/HH5_Save-ArchR]

        // Each stage is now a separate object in the channel,
        // NOTE: stage objects are not reclustered so UMAPs will look weird
        // NOTE: all stage objects should have the same peak set which was calculated on the full data


    } else {
       
       // NEED TO UPDATE
       METADATA_PEAKCALL_PROCESSED( params.peakcall_processed_sample_sheet )
       ch_peakcall_processed = METADATA_PEAKCALL_PROCESSED.out.metadata 

    }

    ///////////////////////////////////////////////////////////////////////////
    /////////////////////    SINGLE CELL PROCESSING      //////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // reclusters and calls peaks on stages using ArchR 
    // integrates stages at single cell level using ArchR

    if(!skip_singlecell_processing){

        // Extract just the stage objects from upstream processing
        ch_upstream_processed
            .filter{ meta, data -> meta.sample_id != 'FullData'}
            .set{ ch_atac_stages }

        // ch_atac_stages.view()
        // [[sample_id:HH5], [Upstream_processing/FILTERING/HH5/rds_files/HH5_Save-ArchR]]
        // [[sample_id:HH6], [Upstream_processing/FILTERING/HH6/rds_files/HH6_Save-ArchR]]
        // [[sample_id:HH7], [Upstream_processing/FILTERING/HH7/rds_files/HH7_Save-ArchR]]
        // [[sample_id:ss4], [Upstream_processing/FILTERING/ss4/rds_files/ss4_Save-ArchR]]
        // [[sample_id:ss8], [Upstream_processing/FILTERING/ss8/rds_files/ss8_Save-ArchR]]

        // re-run clustering and peak calling to get plots for individual stages
        ARCHR_STAGE_PEAKS_WF( ch_atac_stages )

        // read in RNA data
        METADATA_RNA_SC( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                            // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                            // etc
   
        // combine ATAC and RNA data
        ARCHR_STAGE_PEAKS_WF.out.output // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( METADATA_RNA_SC.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
            .groupTuple( by:0 ) // [[sample_id:HH5], [[HH5_Save-ArchR], [HH5_splitstage_data/rds_files/HH5_clustered_data.RDS]]]
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            // .view()
            .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

        // ARCHR: Integrates RNA and ATAC data at single cell level
        ARCHR_INTEGRATING_WF( ch_integrate )  // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]], etc]

        // Run differential accessibility tests between stage-specific peaks
        ARCHR_STAGE_DIFF_PEAKS_WF( ARCHR_INTEGRATING_WF.out )

                /// THINGS TO ADD ///
                // recreate how I found the enhancers I tested

        // Each stage is a separate object in the channel,
        // NOTE: stage objects are reclustered individually
        // NOTE: all stage objects have a different peak set which was calculated on each stage individually
        

    } else {
       
         // NEED TO UPDATE
    //    METADATA_SINGLECELL_PROCESSED( params.singlecell_processed_sample_sheet )
    //    ch_singlecell_processed = METADATA_SINGLECELL_PROCESSED.out.metadata 

    }


    /////////////////////////////////////////////////////////////////////////
    /////////////////////    METACELLS PROCESSING      //////////////////////
    /////////////////////////////////////////////////////////////////////////
    // takes the RNA stage objects and
    // takes the ATAC stage objects which were clustered and peak called altogether (consensus peak set)
    // calculates metacells for RNA and ATAC stages
    // runs integration between RNA and ATAC metacells
    // finds most variable peaks and clusters them into modules

    if(!skip_metacell_processing){

        METADATA_PEAKCALL_PROCESSED( params.peakcall_processed_sample_sheet )
        ch_peakcall_processed = METADATA_PEAKCALL_PROCESSED.out.metadata 
        // ch_peakcall_processed.view()
        // [[sample_id:HH5], [FullData/Split_stages/rds_files/HH5_Save-ArchR]]
        // [[sample_id:HH6], [FullData/Split_stages/rds_files/HH6_Save-ArchR]]
        // [[sample_id:HH7], [FullData/Split_stages/rds_files/HH7_Save-ArchR]]
        // [[sample_id:ss4], [FullData/Split_stages/rds_files/ss4_Save-ArchR]]
        // [[sample_id:ss8], [FullData/Split_stages/rds_files/ss8_Save-ArchR]]

        ///////     Calculate SEACells      ///////

        // Run Metacells on ATAC stages
        SEACELLS_ATAC_WF( ch_peakcall_processed, ch_binary_knowledge_matrix )
             
        //read in RNA data (stages only)
        METADATA_RNA( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                                // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                                // etc
        // Run Metacells on RNA stages
        SEACELLS_RNA_WF( METADATA_RNA.out.metadata, ch_binary_knowledge_matrix )

        ///////     Integrate SEACells      ///////

        // will these different outputs channel in stage by stage??
        SEACELLS_INTEGRATING_WF( SEACELLS_RNA_WF.out.seacells_anndata_processed_classified, SEACELLS_ATAC_WF.out.seacells_anndata_processed_classified, SEACELLS_ATAC_WF.out.seacells_seurat_processed_classified, SEACELLS_ATAC_WF.out.seacell_outputs_named )

        ///////     Cluster peaks      ///////
        CLUSTER_PEAKS_WF( SEACELLS_ATAC_WF.out.seacell_outputs_named, SEACELLS_INTEGRATING_WF.out.processed_integration_output )

        ///////     Transfer SEACells labels onto single cells      ///////
        //SEACELLS_ATAC_WF.out.seacell_outputs_named.view()

        METADATA_SINGLECELL_PROCESSED( params.singlecell_processed_sample_sheet )
        ch_singlecell_processed = METADATA_SINGLECELL_PROCESSED.out.metadata 

        // ch_singlecell_processed.view()
        // [[sample_id:HH5], [HH5/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH5_Save-ArchR]]
        // [[sample_id:HH6], [HH6/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH6_Save-ArchR]]
        // [[sample_id:HH7], [HH7/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH7_Save-ArchR]]
        // [[sample_id:ss4], [ss4/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss4_Save-ArchR]]
        // [[sample_id:ss8], [ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR]]

        // and make Txdb object for plotting - at some point just save the TxDB object saved in the first preprocessing step instead
        // Channel
        //     .value(params.reference)
        //     .map { row -> [[sample_id:'dummy'], row] }
        //     .set{ch_dummy}
        //[[sample_id:dummy], /nemo/lab/briscoej/working/hamrude/raw_data/genomes/galgal6/tag_chroms.gtf]

        // MAKE_TXDB(ch_dummy)

        // EXTRACT_EXONS(ch_gtf)

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////    INTEGRATE CLUSTERS AND METACELLS AT SINGLE CELL LEVEL: STAGES   //////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // takes the metacell ID : single cell ID dictionary
    // takes the consensus peak set GRanges object written as a .csv
    // takes the individually processed single cell data
    // transfers the consensus peakset onto each stage of the ATAC data object 
    // transfer the metacell IDs + integrated labels onto each stage of the ATAC data object

    // runs differential peak analysis on metacells to see how that differs from clusters (using stage-specific peaks)



    
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
