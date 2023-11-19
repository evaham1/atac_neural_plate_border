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

// include { METADATA as METADATA_PEAKCALL_PROCESSED } from "$baseDir/subworkflows/local/metadata"
// include { METADATA as METADATA_METACELL_CSVS } from "$baseDir/subworkflows/local/metadata"
// include { METADATA as METADATA_MEGA_INPUT } from "$baseDir/subworkflows/local/metadata"



// UPSTREAM PROCESSING
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PREPROCESSING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Filtering"

// SINGLE CELL PROCESSING
include { METADATA as METADATA_UPSTREAM_PROCESSED } from "$baseDir/subworkflows/local/metadata"

include {R as CLUSTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_peak_calling.R", checkIfExists: true) )
include {R as INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_constrained_integration_coaccessibility.R", checkIfExists: true) )
include {R as MOTIF_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_motif_analysis.R", checkIfExists: true) )
include {R as CLUSTER_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include { METADATA as METADATA_RNA_SC } from "$baseDir/subworkflows/local/metadata"
include {R as TRANSFER_LABELS_AND_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_transfer_peaks_and_labels.R", checkIfExists: true) )
include {R as MOTIF_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_motif_analysis.R", checkIfExists: true) )

include {R as REMOVE_CONTAM_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_subsetting.R", checkIfExists: true) )
include {R as RECLUSTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include { METADATA as METADATA_RNA_LATENT_TIME } from "$baseDir/subworkflows/local/metadata"
include {R as TRANSFER_LATENT_TIME} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/Transfer_latent_time.R", checkIfExists: true) )
include {R as TRANSFER_LATENT_TIME_MINUS_CONTAM} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/Transfer_latent_time.R", checkIfExists: true) )

include {R as PLOT_DIFF_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_plot_diff_peaks.R", checkIfExists: true) )
include {R as PLOT_DIM_RED_GENOMIC_SUBSETS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_dim_red_genomic_subsets.R", checkIfExists: true) )

include {R as PLOT_MOTIF_CLUSTERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_plot_motifs.R", checkIfExists: true) )
include {R as PLOT_COACCESSIBILITY_CLUSTERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_plot_coaccessibility.R", checkIfExists: true) )

// MEGA PROCESSING
include { METADATA as METADATA_SINGLECELL_PROCESSED } from "$baseDir/subworkflows/local/metadata"

include {R as REMOVE_HH4} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/remove_HH4_RNA_data.R", checkIfExists: true) )
include {R as ARCHR_TO_SEURAT} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/ArchR_to_seurat.R", checkIfExists: true) )
include {R as MEGA_PAIRING_CHROMVAR} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_make_multiome.R", checkIfExists: true) )
include {R as MEGA_GRNI} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_GRNi.R", checkIfExists: true) )
include {R as MEGA_GRN_VIS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_GRN_vis.R", checkIfExists: true) )
include {R as MEGA_GRNI_GMS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_GRNi_GMs.R", checkIfExists: true) )
include {R as MEGA_GRN_GMS_VIS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_GRN_vis.R", checkIfExists: true) )


// METACELL PROCESSING
include { METADATA as METADATA_METACELL_INPUT } from "$baseDir/subworkflows/local/metadata"

include { SEACELLS_ATAC_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_ATAC_WF"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { SEACELLS_RNA_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_RNA_WF"
include { SEACELLS_INTEGRATING_WF } from "$baseDir/subworkflows/local/PROCESSING/SEACells_integration_WF"

//PEAK CLUSTERING
include { CLUSTER_PEAKS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/cluster_peaks_WF"

// // MISC
// include {R as MAKE_TXDB} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/gtf_to_txdb.R", checkIfExists: true) )
// include {EXTRACT_EXONS} from "$baseDir/modules/local/extract_exons/main"

// // DOWNSTREAM PROCESSING WORKFLOWS ~ MULTIVIEW
include {R as MOTIF_ANALYSIS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_motif_analysis.R", checkIfExists: true) )

// include {R as TRANSFER_METACELL_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/ATAC_seacell_purity.R", checkIfExists: true) )
// include {R as TRANSFER_CONSENSUS_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_transfer_peaks.R", checkIfExists: true) )
// include {R as PLOT_DIFF_PEAKS_METACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_peaks/diff_peaks_plots.R", checkIfExists: true) )


// 
// include { FIND_ENHANCERS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/find_enhancers_WF"
//include { COMPARE_VARIABILITY } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_compare_variability"
//include { NPB_SUBSET } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_npb_subset"

// PARAMS FOR PIPELINE SWITCHING
def skip_upstream_processing = params.skip_upstream_processing ? true : false
def skip_sc_processing = params.skip_sc_processing ? true : false
def skip_metacell_processing = params.skip_metacell_processing ? true : false
def skip_multiview_processing = params.skip_multiview_processing ? true : false
def skip_mega_processing = params.skip_mega_processing ? true : false


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

// set channel to just fasta file
Channel
    .value(params.fasta)
    .set{ch_fasta}

// set channel to P2G linkage csv file created using integration wf
Channel
    .value(params.p2g)
    .set{ch_p2g}

// set channel to seurat object
Channel
    .value(params.seurat)
    .set{ch_seurat}

// set channel to atac object to do footprinting for key factors
Channel
    .value(params.atac)
    .set{ch_atac}

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
    /////////////////////      SINGLE CELL PROCESSING      ///////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    // aim of this section is to end up with a FullData object and stage objects that have:
            // clusters at appropriate resolution
            // a consensus peak set (calculated on the full data)
            // consensus cell state labels transferred from the RNA data (calculated on the stages)
            // for full data also have minus contam and plus latent time values

    if(!skip_sc_processing){

        ////    Process full data   ////
        // Extract just the full data object
        ch_upstream_processed
            .filter{ meta, data -> meta.sample_id == 'FullData'}
            .set{ ch_full }

        // Cluster full data
        CLUSTER_FULL( ch_full )

        // Call peaks on full data
        PEAK_CALL( CLUSTER_FULL.out )

        // Read in full RNA data
        METADATA_RNA_SC( params.rna_fulldata_sample_sheet ) //[[sample_id:FullData], [seurat_label_transfer.RDS]]

        // Combine ATAC and RNA full data
        PEAK_CALL.out // [[sample_id:FullData], [/ArchRLogs, /Rplots.pdf, /csv_files, /plots, /rds_files, /tmp]]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( METADATA_RNA_SC.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
            .groupTuple( by:0 ) // [[sample_id:HH5], [[HH5_Save-ArchR], [HH5_splitstage_data/rds_files/HH5_clustered_data.RDS]]]
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            //.view() //[[sample_id:FullData], [rds_files, seurat_label_transfer.RDS]]
            .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

        // Integrate RNA and ATAC data full data
        INTEGRATE( ch_integrate )

        // // Run motif analysis on full data
        // MOTIF_FULL( INTEGRATE.out )

        ////    Process stage data   ////
        // Extract just the stage data objects
        ch_upstream_processed
            .filter{ meta, data -> meta.sample_id != 'FullData'}
            .set{ ch_stages }

        // Cluster individual stages
        CLUSTER_STAGES( ch_stages )

        ////    Transfer labels and peak data from full data onto stages Data   ////
        INTEGRATE.out
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            //.view() FullData_Save-ArchR
            .set{ ch_full }
        CLUSTER_STAGES.out
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            // .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            // .view()
            //     [[sample_id:HH6], [rds_files]]
            //     [[sample_id:HH5], [rds_files]]
            //     [[sample_id:HH7], [rds_files]]
            //     [[sample_id:ss8], [rds_files]]
            //     [[sample_id:ss4], [rds_files]]
            .set{ stages_data }
        stages_data
            .combine(ch_full)
            //.view() //[[sample_id:ss8], [ss8_Save-ArchR], FullData_Save-ArchR]
            .map { row -> [row[0], [row[1][0], row[2]]]}
            //.view() //[[sample_id:ss8], [ss8_Save-ArchR, FullData_Save-ArchR]]
            .set{ch_transfer}
        TRANSFER_LABELS_AND_PEAKS(ch_transfer)

        // // Run motif analysis on stages data
        // MOTIF_STAGES( TRANSFER_LABELS_AND_PEAKS.out )

        //    Extra processing with full data  ////
        
        // Read in RNA object with latent time
        METADATA_RNA_LATENT_TIME( params.rna_latent_time_sample_sheet )

        // Transfer latent time from RNA full data to ATAC full data
        INTEGRATE.out
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( METADATA_RNA_LATENT_TIME.out.metadata )
            .groupTuple( by:0 )
            //.view() //[[sample_id:FullData], [[rds_files], [seurat_label_transfer_latent_time.RDS]]]
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            //.view() //[[sample_id:FullData], [rds_files, seurat_label_transfer_latent_time.RDS]]
            .set {ch_transfer_latent_time} //[[sample_id:FullData], [plots, rds_files]]
        TRANSFER_LATENT_TIME( ch_transfer_latent_time )

        // Remove contam from Full data and re-cluster
        REMOVE_CONTAM_FULL( INTEGRATE.out )
        RECLUSTER_FULL( REMOVE_CONTAM_FULL.out )

        // Transfer latent time from RNA full data to ATAC full data WITHOUT CONTAM
        RECLUSTER_FULL.out
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( METADATA_RNA_LATENT_TIME.out.metadata )
            .groupTuple( by:0 )
            //.view() //[[sample_id:FullData], [[rds_files], [seurat_label_transfer_latent_time.RDS]]]
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            //.view() //[[sample_id:FullData], [rds_files, seurat_label_transfer_latent_time.RDS]]
            .set {ch_transfer_latent_time_minus_contam} //[[sample_id:FullData], [plots, rds_files]]
        TRANSFER_LATENT_TIME_MINUS_CONTAM( ch_transfer_latent_time_minus_contam )


        ////// PLOTTING /////// - maybe move to multiview section once I've tested that it works

        // Plot diff peaks between clusters in stages
        PLOT_DIFF_PEAKS( TRANSFER_LABELS_AND_PEAKS.out )

        // Plot UMAPs of dim reduction using different subsets
        PLOT_DIM_RED_GENOMIC_SUBSETS( TRANSFER_LABELS_AND_PEAKS.out )

        // Coaccessibility plots grouped by clusters
        // need to have the stages objects + the coaccessibility csv from full data
        INTEGRATE.out
            .map{it[1].findAll{it =~ /csv_files/}[0]}
            .set{ ch_full_coaccessibility }
        TRANSFER_LABELS_AND_PEAKS.out
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .set{ stages_labelled_data }
        stages_labelled_data
            .combine(ch_full_coaccessibility)
            //.view() //[[sample_id:ss8], [rds_files], csv_files]
            .map { row -> [row[0], [row[1][0], row[2]]]}
            //.view() //[[sample_id:HH7], [rds_files, csv_files]]
            .set{ch_stages_coaccessibility}
        PLOT_COACCESSIBILITY_CLUSTERS( ch_stages_coaccessibility )
        
        // Motif analysis plots grouped by clusters
        // just needs the archr objects
        // PLOT_MOTIF_CLUSTERS( MOTIF_STAGES.out )
        




    } else {
       
       // NEED TO UPDATE
    //    METADATA_SC_PROCESSED( params.peakcall_processed_sample_sheet )
    //    ch_peakcall_processed = METADATA_PEAKCALL_PROCESSED.out.metadata 

    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////    MEGA PROCESSING      //////////////////////
    /////////////////////////////////////////////////////////////////////////

    if(!skip_mega_processing){

        ////    Prep ATAC data    ////

        // read in ATAC full transfer label object with latent time
        METADATA_SINGLECELL_PROCESSED( params.mega_input_sample_sheet ) // single cell data with individual peaks called
        ch_singlecell_processed = METADATA_SINGLECELL_PROCESSED.out.metadata
        ch_singlecell_processed.view()

        // convert ArchR full data object into seurat object
        ARCHR_TO_SEURAT( ch_singlecell_processed )

        ////    Prep RNA data    ////

        // read in RNA full data
        METADATA_RNA_SC( params.rna_latent_time_sample_sheet )
        //remove HH4 from RNA full data
        REMOVE_HH4( METADATA_RNA_SC.out.metadata )
        // extract RNA seurat object
        REMOVE_HH4.out //[[sample_id:FullData], [plots, rds_files]]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .flatMap {it[1][0].listFiles()}
            //.view() //seurat_label_transfer_minus_HH4.RDS
            .map { row -> [[sample_id:'FullData'], row] }
            .set { ch_rna }
   
        ////    Integrate RNA and ATAC data    ////

        // combine ATAC and RNA channels
        ARCHR_TO_SEURAT.out // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            //.view() //[[sample_id:FullData], [rds_files]]
            .concat( ch_rna )
            .groupTuple( by:0 )
            //.view() //[ [sample_id:FullData], [[rds_files], seurat_label_transfer_minus_HH4.RDS] ]
            .map{ [ it[0], [ it[1][0][0], it[1][1] ] ] }
            //.view() //[[sample_id:FullData], [rds_files, seurat_label_transfer_minus_HH4.RDS]]
            .set {ch_integrate} //[[sample_id:FullData], [plots, rds_files]]

        // use previously calculated integration to pair ATAC and RNA cells -> fake multimodal data
        MEGA_PAIRING_CHROMVAR( ch_integrate )

        // then run scMEGA GRNi on the full data paired seurat object
        MEGA_PAIRING_CHROMVAR.out
            .combine(ch_p2g)
            .map{[it[0], it[1] + it[2]]}
            //.view()
            .set {ch_grni} // add p2g csv file
        MEGA_GRNI( ch_grni )
        MEGA_GRNI.out
            .combine(ch_seurat)
            .combine(ch_atac)
            .map{[it[0], it[1] + it[2] + it[3]]}
            .view()
            .set {ch_grni_vis} // add p2g csv file
        MEGA_GRN_VIS( ch_grni_vis )

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

        // read in stages objects with consensus peak set
        METADATA_METACELL_INPUT( params.metacell_input_sample_sheet )
        ch_metacell_input = METADATA_METACELL_INPUT.out.metadata 
        //ch_metacell_input.view()
        // [[sample_id:HH5], [FullData/Split_stages/rds_files/HH5_Save-ArchR]]
        // [[sample_id:HH6], [FullData/Split_stages/rds_files/HH6_Save-ArchR]]
        // [[sample_id:HH7], [FullData/Split_stages/rds_files/HH7_Save-ArchR]]
        // [[sample_id:ss4], [FullData/Split_stages/rds_files/ss4_Save-ArchR]]
        // [[sample_id:ss8], [FullData/Split_stages/rds_files/ss8_Save-ArchR]]

        ///////     Calculate SEACells      ///////

        // Run Metacells on ATAC stages
        SEACELLS_ATAC_WF( ch_metacell_input, ch_binary_knowledge_matrix )
             
        //read in RNA data (stages only)
        METADATA_RNA( params.rna_stages_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                                // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                                // etc
        // Run Metacells on RNA stages
        SEACELLS_RNA_WF( METADATA_RNA.out.metadata, ch_binary_knowledge_matrix )
        // maybe in here repurpose the metacell purity script to check they are pure in the RNA data?

        ///////     Integrate SEACells      ///////

        // will these different outputs channel in stage by stage??
        SEACELLS_INTEGRATING_WF( SEACELLS_RNA_WF.out.seacells_anndata_processed_classified, SEACELLS_ATAC_WF.out.seacells_anndata_processed_classified, SEACELLS_ATAC_WF.out.seacells_seurat_processed, SEACELLS_ATAC_WF.out.seacell_outputs_named )
        // maybe in here add schelper cell type broad labels

        ///////     Cluster peaks      ///////
        CLUSTER_PEAKS_WF( SEACELLS_ATAC_WF.out.seacell_outputs_named, SEACELLS_INTEGRATING_WF.out.processed_integration_output, ch_fasta )

        // ch_singlecell_processed.view()
        // [[sample_id:HH5], [HH5/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH5_Save-ArchR]]
        // [[sample_id:HH6], [HH6/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH6_Save-ArchR]]
        // [[sample_id:HH7], [HH7/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH7_Save-ArchR]]
        // [[sample_id:ss4], [ss4/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss4_Save-ArchR]]
        // [[sample_id:ss8], [ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR]]

        // run script to transfer metacell IDs to single cells on each ArchR stage object - script made 'ArchR_seacell_purity'

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////    INTEGRATE CLUSTERS AND METACELLS AT SINGLE CELL LEVEL: STAGES   //////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // takes the metacell ID : single cell ID dictionary
    // takes the consensus peak set GRanges object written as a .csv
    // takes the individually processed single cell data
    // transfers the consensus peakset onto each stage of the ATAC data object 
    // transfer the metacell IDs + integrated labels onto each stage of the ATAC data object

    if(!skip_multiview_processing){

        // METADATA_PEAKCALL_PROCESSED( params.peakcall_processed_sample_sheet ) // single cell data with consensus peaks called
        // ch_peakcall_processed = METADATA_PEAKCALL_PROCESSED.out.metadata 

        METADATA_SINGLECELL_PROCESSED( params.singlecell_processed_sample_sheet ) // single cell data with individual peaks called
        ch_singlecell_processed = METADATA_SINGLECELL_PROCESSED.out.metadata

        // METADATA_METACELL_CSVS( params.metacell_csvs_sample_sheet ) // csv files with metacell IDs
        // ch_metadata_csvs = METADATA_METACELL_CSVS.out.metadata
        // ch_metadata_csvs.view()
        // [[sample_id:HH5], [HH5/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/HH5_cell_metadata.csv]]
        // [[sample_id:HH6], [HH6/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/HH6_cell_metadata.csv]]
        // [[sample_id:HH7], [HH7/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/HH7_cell_metadata.csv]]
        // [[sample_id:ss4], [ss4/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/ss4_cell_metadata.csv]]
        // [[sample_id:ss8], [ss8/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/ss8_cell_metadata.csv]]


        ///////     Transfer SEACells labels onto single cells      ///////
        // and check how they correspond with other single cell labels - script made 'ArchR_seacell_purity'
        // run peak calling and diff peak analysis to see how this compares to cluster-level analysis (just rerun ARCHR_STAGE_DIFF_PEAKS_WF)

        // // combine ArchR objects and metacell csvs
        // ch_singlecell_processed
        //     .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
        //     .concat( ch_metadata_csvs )
        //     .groupTuple( by:0 )
        //     .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
        //     //.view()
        //     .set {ch_transfer_metacell_IDs}
        // // [[sample_id:HH5], [/HH5/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH5_Save-ArchR, /HH5/Integrated_SEACells_label_transfer/rds_files/HH5_ATAC_singlecell_integration_map.csv]]
        // // [[sample_id:HH6], [/HH6/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH6_Save-ArchR, /HH6/Integrated_SEACells_label_transfer/rds_files/HH6_ATAC_singlecell_integration_map.csv]]
        // // [[sample_id:HH7], [/HH7/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH7_Save-ArchR, /HH7/Integrated_SEACells_label_transfer/rds_files/HH7_ATAC_singlecell_integration_map.csv]]
        // // [[sample_id:ss4], [/ss4/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss4_Save-ArchR, /ss4/Integrated_SEACells_label_transfer/rds_files/ss4_ATAC_singlecell_integration_map.csv]]
        // // [[sample_id:ss8], [/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR, /ss8/Integrated_SEACells_label_transfer/rds_files/ss8_ATAC_singlecell_integration_map.csv]]
                
        // // run script to transfer metacell IDs to single cells on each ArchR stage object - script made 'ArchR_seacell_purity'
        // TRANSFER_METACELL_LABELS( ch_transfer_metacell_IDs )

        // // call peaks on the metacell integrated labels and visualise their differential accessibility (to comporate to cluster analysis)
        // //PEAK_CALL_METACELLS( TRANSFER_METACELL_LABELS.out )
        // PLOT_DIFF_PEAKS_METACELLS( PEAK_CALL_METACELLS.out )


        // Take peaks from PMs and run them through Alex's pipeline??


        ///////     Visualisations on final data      ///////
        // Motif analysis

        MOTIF_ANALYSIS( ch_singlecell_processed )


    }

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
