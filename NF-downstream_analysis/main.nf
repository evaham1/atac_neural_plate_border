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
include { METADATA as METADATA_METACELL_CSVS } from "$baseDir/subworkflows/local/metadata"
include { METADATA as METADATA_MEGA_INPUT } from "$baseDir/subworkflows/local/metadata"

// UPSTREAM PROCESSING
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PREPROCESSING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Filtering"

// PROCESSING - FULL DATA
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_peak_calling.R", checkIfExists: true) )
include {R as CO_ACCESSIBILITY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_coaccessibility.R", checkIfExists: true) )
include {R as SPLIT_STAGES_PROCESSED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_split_stages.R", checkIfExists: true) )

// PROCESSING - STAGES
include { METADATA as METADATA_RNA_SC } from "$baseDir/subworkflows/local/metadata"
include { ARCHR_INTEGRATING_WF } from "$baseDir/subworkflows/local/PROCESSING/archr_integration_WF"
include {R as PLOT_DIFF_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_peaks/diff_peaks_plots.R", checkIfExists: true) )
include {R as MOTIF_ANALYSIS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_motif_analysis.R", checkIfExists: true) )

// MEGA PROCESSING
include {R as ARCHR_TO_SEURAT} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/ArchR_to_seurat.R", checkIfExists: true) )
include {R as MEGA_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_integration.R", checkIfExists: true) )
include {R as MEGA_PAIRING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_pairing.R", checkIfExists: true) )
include {R as MEGA_GRNI} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/scMEGA/MEGA_GRNi.R", checkIfExists: true) )
include {R as REMOVE_HH4} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/remove_HH4_RNA_data.R", checkIfExists: true) )

// METACELL PROCESSING
include { SEACELLS_ATAC_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_ATAC_WF"
include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { SEACELLS_RNA_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_RNA_WF"
include { SEACELLS_INTEGRATING_WF } from "$baseDir/subworkflows/local/PROCESSING/SEACells_integration_WF"

//PEAK CLUSTERING
include { CLUSTER_PEAKS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/cluster_peaks_WF"

// // MISC
// include {R as MAKE_TXDB} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/data_conversion/gtf_to_txdb.R", checkIfExists: true) )
// include {EXTRACT_EXONS} from "$baseDir/modules/local/extract_exons/main"

// DOWNSTREAM PROCESSING WORKFLOWS ~ MULTIVIEW
include {R as TRANSFER_METACELL_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/ATAC_seacell_purity.R", checkIfExists: true) )
include {R as TRANSFER_CONSENSUS_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_transfer_peaks.R", checkIfExists: true) )
include {R as PLOT_DIFF_PEAKS_METACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_peaks/diff_peaks_plots.R", checkIfExists: true) )


// 
// include { FIND_ENHANCERS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/find_enhancers_WF"
//include { COMPARE_VARIABILITY } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_compare_variability"
//include { NPB_SUBSET } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_npb_subset"

// PARAMS FOR PIPELINE SWITCHING
def skip_upstream_processing = params.skip_upstream_processing ? true : false
def skip_peakcall_processing = params.skip_peakcall_processing ? true : false
def skip_metacell_processing = params.skip_metacell_processing ? true : false
def skip_singlecell_processing = params.skip_singlecell_processing ? true : false
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

        // Cluster, call peaks and calculate interactions by accessibility on full data
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
        METADATA_PEAKCALL_PROCESSED( params.peakcall_processed_sample_sheet )
        ch_atac_stages = METADATA_PEAKCALL_PROCESSED.out.metadata 
        // [[sample_id:HH5], [FullData/Split_stages/rds_files/HH5_Save-ArchR]]
        // [[sample_id:HH6], [FullData/Split_stages/rds_files/HH6_Save-ArchR]]
        // [[sample_id:HH7], [FullData/Split_stages/rds_files/HH7_Save-ArchR]]
        // [[sample_id:ss4], [FullData/Split_stages/rds_files/ss4_Save-ArchR]]
        // [[sample_id:ss8], [FullData/Split_stages/rds_files/ss8_Save-ArchR]]

        // re-run clustering - keep peaks from full data (double check this works ok when running differential peaks)
        CLUSTER( ch_atac_stages )

        // read in RNA data
        METADATA_RNA_SC( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                            // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                            // etc
   
        // combine ATAC and RNA data
        CLUSTER.out // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( METADATA_RNA_SC.out.metadata ) // [ [sample_id:HH5], [HH5_clustered_data.RDS] ]
            .groupTuple( by:0 ) // [[sample_id:HH5], [[HH5_Save-ArchR], [HH5_splitstage_data/rds_files/HH5_clustered_data.RDS]]]
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            // .view()
            .set {ch_integrate} //[ [sample_id:HH5], [HH5_Save-ArchR, HH5_clustered_data.RDS] ]

        // ARCHR: Integrates RNA and ATAC data at single cell level
        ARCHR_INTEGRATING_WF( ch_integrate )  // [ [[meta: HH5], [RNA, ATAC]] , [[meta: HH6], [RNA, ATAC]], etc]
            // IN HERE NEED TO EDIT SO SCHELPER_CELL_TYPE_OLD = SCHELPER_CELL_TYPE
            // ALSO MAYBE SHOULD ADD SCHELPER_CELL_TYPE_BROAD LABELS
                    // check these changes have worked!

        // Run motif analysis and footprinting for key TFs
        //MOTIF_ANALYSIS( ARCHR_INTEGRATING_WF.out )
        //CO_ACCESSIBILITY( ARCHR_INTEGRATING_WF.out )

        // Run differential accessibility tests between consensus peaks - check that this is still ok
        //PLOT_DIFF_PEAKS( ARCHR_INTEGRATING_WF.out )
                /// THINGS TO ADD INTO THIS WORKFLOW ///
                // recreate how I found the enhancers I tested (may need to re-call peaks on stages)

        // Each stage is a separate object in the channel,
        // NOTE: stage objects are reclustered individually
        // NOTE: all stage objects have the same peak set (consensus full data peakset)
        

    } else {
       
         // NEED TO UPDATE
    //    METADATA_SINGLECELL_PROCESSED( params.singlecell_processed_sample_sheet )
    //    ch_singlecell_processed = METADATA_SINGLECELL_PROCESSED.out.metadata 

    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////    MEGA PROCESSING      //////////////////////
    /////////////////////////////////////////////////////////////////////////

    if(!skip_mega_processing){

        METADATA_SINGLECELL_PROCESSED( params.mega_input_sample_sheet ) // single cell data with individual peaks called
        ch_singlecell_processed = METADATA_SINGLECELL_PROCESSED.out.metadata 

        // ch_singlecell_processed.view()
        //[[sample_id:FullData], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/output/NF-downstream_analysis/Processing/FullData/TransferLabels/rds_files/TransferLabels_Save-ArchR]]

        // re-run clustering - keep peaks from full data (double check this works ok when running differential peaks)
        //CLUSTER( ch_atac_stages )

        // convert ArchR objects into seurat objects
        ARCHR_TO_SEURAT( ch_singlecell_processed )
        // ARCHR_TO_SEURAT.out.view()
        // [[sample_id:HH6], [ArchRLogs, plots, rds_files]]
        // [[sample_id:HH5], [ArchRLogs, plots, rds_files]]
        // [[sample_id:HH7], [ArchRLogs, plots, rds_files]]
        // [[sample_id:ss4], [ArchRLogs, plots, rds_files]]
        // [[sample_id:ss8], [ArchRLogs, plots, rds_files]]

        // read in RNA data
        METADATA_RNA_SC( params.rna_sample_sheet ) // [[sample_id:HH5], [HH5_clustered_data.RDS]]
                                            // [[sample_id:HH6], [HH6_clustered_data.RDS]]
                                            // etc
        // remove HH4 from RNA data
        REMOVE_HH4( METADATA_RNA_SC.out.metadata )
        
        // extract RNA seurat object
        REMOVE_HH4.out //[[sample_id:FullData], [plots, rds_files]]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .flatMap {it[1][0].listFiles()}
            //.view() //seurat_label_transfer_minus_HH4.RDS
            .map { row -> [[sample_id:'FullData'], row] }
            .set { ch_rna }
   
        // combine ATAC and RNA data
        ARCHR_TO_SEURAT.out // [ [sample_id:HH5], [ArchRLogs, Rplots.pdf, plots, rds_files] ]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            //.view() //[[sample_id:FullData], [rds_files]]
            .concat( ch_rna )
            .groupTuple( by:0 )
            //.view() //[ [sample_id:FullData], [[rds_files], seurat_label_transfer_minus_HH4.RDS] ]
            .map{ [ it[0], [ it[1][0][0], it[1][1] ] ] }
            //.view() //[[sample_id:FullData], [rds_files, seurat_label_transfer_minus_HH4.RDS]]
            .set {ch_integrate} //[[sample_id:FullData], [plots, rds_files]]

        // integrate the stages into a coembedding seurat object
        MEGA_INTEGRATION( ch_integrate )

        // use previously calculated integration to pair ATAC and RNA cells -> fake multimodal data
        MEGA_PAIRING( MEGA_INTEGRATION.out )
        
        // then run scMEGA GRNi on the full data paired seurat object
        //MEGA_GRNI( MEGA_INTEGRATION.out )
        

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
        // maybe in here repurpose the metacell purity script to check they are pure in the RNA data?

        ///////     Integrate SEACells      ///////

        // will these different outputs channel in stage by stage??
        SEACELLS_INTEGRATING_WF( SEACELLS_RNA_WF.out.seacells_anndata_processed_classified, SEACELLS_ATAC_WF.out.seacells_anndata_processed_classified, SEACELLS_ATAC_WF.out.seacells_seurat_processed_classified, SEACELLS_ATAC_WF.out.seacell_outputs_named )
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

        METADATA_METACELL_CSVS( params.metacell_csvs_sample_sheet ) // csv files with metacell IDs
        ch_metadata_csvs = METADATA_METACELL_CSVS.out.metadata
        ch_metadata_csvs.view()
        // [[sample_id:HH5], [HH5/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/HH5_cell_metadata.csv]]
        // [[sample_id:HH6], [HH6/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/HH6_cell_metadata.csv]]
        // [[sample_id:HH7], [HH7/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/HH7_cell_metadata.csv]]
        // [[sample_id:ss4], [ss4/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/ss4_cell_metadata.csv]]
        // [[sample_id:ss8], [ss8/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/ss8_cell_metadata.csv]]


        ///////     Transfer SEACells labels onto single cells      ///////
        // and check how they correspond with other single cell labels - script made 'ArchR_seacell_purity'
        // run peak calling and diff peak analysis to see how this compares to cluster-level analysis (just rerun ARCHR_STAGE_DIFF_PEAKS_WF)

        // combine ArchR objects and metacell csvs
        ch_singlecell_processed
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .concat( ch_metadata_csvs )
            .groupTuple( by:0 )
            .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
            //.view()
            .set {ch_transfer_metacell_IDs}
        // [[sample_id:HH5], [/HH5/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH5_Save-ArchR, /HH5/Integrated_SEACells_label_transfer/rds_files/HH5_ATAC_singlecell_integration_map.csv]]
        // [[sample_id:HH6], [/HH6/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH6_Save-ArchR, /HH6/Integrated_SEACells_label_transfer/rds_files/HH6_ATAC_singlecell_integration_map.csv]]
        // [[sample_id:HH7], [/HH7/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH7_Save-ArchR, /HH7/Integrated_SEACells_label_transfer/rds_files/HH7_ATAC_singlecell_integration_map.csv]]
        // [[sample_id:ss4], [/ss4/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss4_Save-ArchR, /ss4/Integrated_SEACells_label_transfer/rds_files/ss4_ATAC_singlecell_integration_map.csv]]
        // [[sample_id:ss8], [/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR, /ss8/Integrated_SEACells_label_transfer/rds_files/ss8_ATAC_singlecell_integration_map.csv]]
                
        // run script to transfer metacell IDs to single cells on each ArchR stage object - script made 'ArchR_seacell_purity'
        TRANSFER_METACELL_LABELS( ch_transfer_metacell_IDs )

        // call peaks on the metacell integrated labels and visualise their differential accessibility (to comporate to cluster analysis)
        //PEAK_CALL_METACELLS( TRANSFER_METACELL_LABELS.out )
        PLOT_DIFF_PEAKS_METACELLS( PEAK_CALL_METACELLS.out )

        ///////     Transfer full data peak set onto individual stages      ///////
        // shouldn't have to do this if always do use the consensus peak set
        // combine ArchR object with metacell IDs and consensus peak set
        // TRANSFER_METACELL_LABELS.out
        //     .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
        //     .concat( ch_peakcall_processed )
        //     .groupTuple( by:0 )
        //     .map{ [ it[0], [ it[1][0][0], it[1][1][0] ] ] }
        //     .view()
        //     .set {ch_transfer_peaks}

        // // makes it easier to work with as can merge and split stages at will
        // TRANSFER_CONSENSUS_PEAKS( ch_transfer_peaks )


        // Take peaks from PMs and run them through Alex's pipeline??

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
