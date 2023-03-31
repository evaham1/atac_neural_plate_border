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

// UPSTREAM PROCESSING WORKFLOWS
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include { PREPROCESSING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Preprocessing"
include { FILTERING } from "$baseDir/subworkflows/local/UPSTREAM_PROCESSING/Filtering"

// PROCESSING WORKFLOWS AND MODULES
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as SPLIT_STAGES_PROCESSED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_split_stages.R", checkIfExists: true) )

//CALCULATING SEACELL WFs
include { SEACELLS_ATAC_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_ATAC_WF"

include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { SEACELLS_RNA_WF } from "$baseDir/subworkflows/local/PROCESSING/seacells_RNA_WF"

//INTEGRATING SEACELLS
include {PYTHON as INTEGRATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/Integration/SEACells_integration.py", checkIfExists: true) )

//PEAK CLUSTERING
include {R as COMBINE_METACELL_COUNTS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Metacell_processes/Combine_summarised_counts.R", checkIfExists: true) )
include { CLUSTER_PEAKS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/cluster_peaks_WF"


// DOWNSTREAM PROCESSING WORKFLOWS

// 
// include { FIND_ENHANCERS_WF } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/find_enhancers_WF"
//include { COMPARE_VARIABILITY } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_compare_variability"
//include { NPB_SUBSET } from "$baseDir/subworkflows/local/DOWNSTREAM_PROCESSING/archr_npb_subset"

// PARAMS
def skip_upstream_processing = params.skip_upstream_processing ? true : false
def skip_peakcall_processing = params.skip_peakcall_processing ? true : false
def skip_metacell_processing = params.skip_metacell_processing ? true : false



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


    /////////////////////////////////////////////////////////////////
    /////////////////////    PEAK CALLING      //////////////////////
    /////////////////////////////////////////////////////////////////
    // calls peaks on full data object to obtain consensus peak set and peak counts

    if(!skip_peakcall_processing){

        // Extract just the full data object
        ch_upstream_processed
            .filter{ meta, data -> meta.sample_id == 'FullData'}
            .set{ ch_full }

        // Cluster and call peaks on full data
        //ch_full.view() //[[sample_id:FullData], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/output/NF-downstream_analysis/Upstream_processing/FILTERING/FullData/rds_files/FullData_Save-ArchR]]
        CLUSTER( ch_full )
        PEAK_CALL( CLUSTER.out )

        // Split full data into stages
        SPLIT_STAGES_PROCESSED( PEAK_CALL.out )
        SPLIT_STAGES_PROCESSED.out //[[meta], [plots, rds_files]]
            .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
            .flatMap {it[1][0].listFiles()}
            .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
            .set { ch_peakcall_processed }

        // ch_peakcall_processed.view()
        // [[sample_id:HH6], rds_files/HH6_Save-ArchR]
        // [[sample_id:HH7], rds_files/HH7_Save-ArchR]
        // [[sample_id:ss4], rds_files/ss4_Save-ArchR]
        // [[sample_id:ss8], rds_files/ss8_Save-ArchR]
        // [[sample_id:HH5], rds_files/HH5_Save-ArchR]

    } else {
       
       METADATA_PEAKCALL_PROCESSED( params.peakcall_processed_sample_sheet )
       ch_peakcall_processed = METADATA_PEAKCALL_PROCESSED.out.metadata 

    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////    METACELLS PROCESSING      //////////////////////
    /////////////////////////////////////////////////////////////////////////
    // calculates metacells for RNA and ATAC stages
    // run integration between RNA and ATAC metacells

    if(!skip_metacell_processing){

        // Run Metacells on ATAC stages
        SEACELLS_ATAC_WF( ch_peakcall_processed, ch_binary_knowledge_matrix )
             
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

        //ch_RNA.view()
        // [[sample_id:HH6], c8/0355d43ba613bc2159e10758c26572/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:HH5], 61/920e02828bc1d30d01745280da08cb/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:HH7], d0/6b148d2543b6b195b8ee2088fad0b7/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:ss4], fa/1d64469d41d5a69f74d4d706bf5681/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:ss8], 64/bc0c8d3f045086356bc00a7ebde6ff/rds_files/AnnData_RNA.h5ad]
        
        SEACELLS_ATAC_WF.out.seacells_anndata_processed_classified
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()[0]] }
            .set {ch_ATAC}
        // ch_ATAC.view()
        // [[sample_id:HH6], 6d/4467b6b5de1ef8f33c2edc45d871b6/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:HH7], 4f/2aef89dec2263316c6a5645b259204/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:HH5], 78/c6ea5ab5e6e2cc3e0b5ccd6763f948/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:ss4], a4/9d8ea368614ff7badd70b2135e7046/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:ss8], 4a/a274c2d504b915e6a9fb8369584796/rds_files/AnnData_ATAC.h5ad]

        ch_RNA
            .concat( ch_ATAC )
            .groupTuple( by:0 )
            //.view() //[[sample_id:ss4], [fa/1d64469d41d5a69f74d4d706bf5681/rds_files/AnnData_RNA.h5ad, 4a/72711280da408db03e302a5926456a/rds_files/AnnData_ATAC.h5ad]]
            .set {ch_seacells_to_integrate}
        INTEGRATE_SEACELLS( ch_seacells_to_integrate )

        ///////     Check new scHelper_cell_type labels at single cell level on ATAC data      ///////

        // ch_labels_combined = SEACELLS_ATAC_WF.out.seacells_anndata //seacell ATAC ID to single cell ATAC ID map
        //     .concat( INTEGRATE_SEACELLS.out ) //seacell RNA ID to seacell ATAC ID map + transferred scHelper_cell_type label
        //     .concat( ch_peakcall_processed ) //original ArchR ATAC single cell object
        //     .groupTuple( by:0 ) // all 3 outputs need to be grouped by stage
        //     //.view()
        //SEACELL_LABELS_ON_ATAC( ch_labels_combined ) //take all this info and output ArchR ATAC stage object with new labels generated from single cell integration


        ///////     Run peak clustering on full data      ///////

        // SEACELLS_ATAC_WF.out.seacell_outputs_named.view()
        // [[sample_id:HH7], 58/b3df47b5a798e1acaa3666df878bdb/csv_files]
        // [[sample_id:HH6], 41/8c7ec3c7a7c3bb79dcd52e4bad02b9/csv_files]
        // [[sample_id:ss8], a8/f7a307efa1093759a27fd61ea09350/csv_files]
        // [[sample_id:ss4], 71/b91188feb85928e3c2096811627513/csv_files]
        // [[sample_id:HH5], b5/52af06cc3559a5d45a976102bc509f/csv_files]

        // take the exported_data outputs from SEACell_computation of the ATAC
        ch_metacells_combined = SEACELLS_ATAC_WF.out.seacell_outputs_named // Collect csv files from all stages
            //.map{it[1].findAll{it =~ /csv_files/}[0].listFiles()[0]}
            .map{ meta, data -> [data.findAll{it =~ /csv_files/}[0].listFiles()[0]] }
            //.map{ meta, data -> [data.findAll{it =~ /csv_files/}[0]] } //[[sample_id:FullData], [csv_files, csv_files, csv_files, csv_files, csv_files]]
            //.map{it[1].listFiles()}
            .collect()
            .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]

        ch_metacells_combined.view()
        //[[sample_id:FullData], [[HH6_feature_metadata.csv, HH6_summarised_by_metacells_counts.csv, HH6_cell_metadata.csv], [HH7_feature_metadata.csv, HH7_summarised_by_metacells_counts.csv, HH7_cell_metadata.csv], [ss8_cell_metadata.csv, ss8_summarised_by_metacells_counts.csv, ss8_feature_metadata.csv], [ss4_cell_metadata.csv, ss4_summarised_by_metacells_counts.csv, ss4_feature_metadata.csv], [HH5_feature_metadata.csv, HH5_cell_metadata.csv, HH5_summarised_by_metacells_counts.csv]]]
        //COMBINE_METACELL_COUNTS( ch_metacells_combined ) //combine all the summarised counts into one summarised counts file, check all feature metadata the same and write, combine all cell metadata too?
        
        // run peak clustering wf
        //CLUSTER_PEAKS_WF( COMBINE_METACELL_COUNTS.out )

    }








    
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
