#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R script to combine metacell outputs
include {R as COMBINE_METACELL_COUNTS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/Combine_summarised_counts.R", checkIfExists: true) )

// R scripts to filter peaks and then cluster them using Antler
include {R as FILTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/cluster_peaks/1_filter_peaks.R", checkIfExists: true) )
include {R as CLUSTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/cluster_peaks/2_antler_calculate_peak_modules.R", checkIfExists: true) )

// Module to run homer motif enrichment on peaks
include {HOMER_MOTIF_ENRICHMENT} from "$baseDir/modules/local/homer_motif_enrichment/main"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// cluster peaks

workflow CLUSTER_PEAKS_WF {
    take:
    seacell_output //should be SEACELLS_ATAC_WF.out.seacell_outputs_named
    integration_output //SEACELLS_INTEGRATING.out.processed_integration_output
    ch_fasta // path to reference fasta sequence

    main:

    //seacell_output.view()
        // [[sample_id:HH7], 58/b3df47b5a798e1acaa3666df878bdb/csv_files]
        // [[sample_id:HH6], 41/8c7ec3c7a7c3bb79dcd52e4bad02b9/csv_files]
        // [[sample_id:ss8], a8/f7a307efa1093759a27fd61ea09350/csv_files]
        // [[sample_id:ss4], 71/b91188feb85928e3c2096811627513/csv_files]
        // [[sample_id:HH5], b5/52af06cc3559a5d45a976102bc509f/csv_files]

    // take the exported_data outputs from SEACell_computation of the ATAC
    ch_metacells_combined = seacell_output // Collect csv files from all stages
        .map{ meta, data -> data.listFiles() } //list files inside only folder in output (csv_files)
            // [HH6_feature_metadata.csv, HH6_summarised_by_metacells_counts.csv, HH6_cell_metadata.csv]
            // [ss4_cell_metadata.csv, ss4_summarised_by_metacells_counts.csv, ss4_feature_metadata.csv]
            // [ss8_cell_metadata.csv, ss8_summarised_by_metacells_counts.csv, ss8_feature_metadata.csv]
            // [HH7_feature_metadata.csv, HH7_summarised_by_metacells_counts.csv, HH7_cell_metadata.csv]
            // [HH5_feature_metadata.csv, HH5_cell_metadata.csv, HH5_summarised_by_metacells_counts.csv]
        .flatten() //removes square brackets from each array
        .collect() // puts all arrays together
        .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    // integration_output.view()
    // [[sample_id:ss8], [plots, rds_files]]
    // [[sample_id:ss4], [plots, rds_files]]
    // [[sample_id:HH6], [plots, rds_files]]

    ch_integration_combined = integration_output
        .map{ meta, data -> [data.findAll{it =~ /rds_files/}[0].listFiles()] }
        .flatten() //removes square brackets from each array
        .collect() // puts all arrays together
        .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    //ch_integration_combined.view()
    //[[sample_id:FullData], [HH6_filtered_SEACells_integration_map.csv, HH6_ATAC_singlecell_integration_map.csv, HH6_seacells_seurat_integrated.RDS, ss8_filtered_SEACells_integration_map.csv, ss8_ATAC_singlecell_integration_map.csv, ss8_seacells_seurat_integrated.RDS, HH5_ATAC_singlecell_integration_map.csv, HH5_filtered_SEACells_integration_map.csv, HH5_seacells_seurat_integrated.RDS, HH7_ATAC_singlecell_integration_map.csv, HH7_filtered_SEACells_integration_map.csv, HH7_seacells_seurat_integrated.RDS, ss4_seacells_seurat_integrated.RDS, ss4_filtered_SEACells_integration_map.csv, ss4_ATAC_singlecell_integration_map.csv]]

    ch_metacells_combined
        .concat( ch_integration_combined )
        .groupTuple( by:0 )
        .map{ meta, data -> [meta, data.flatten()] }
        .set { ch_combined_input }

    //ch_combined_input.view()
    //[[sample_id:FullData], [HH6_feature_metadata.csv, HH6_summarised_by_metacells_counts.csv, HH6_cell_metadata.csv, HH7_feature_metadata.csv, HH7_summarised_by_metacells_counts.csv, HH7_cell_metadata.csv, ss4_cell_metadata.csv, ss4_summarised_by_metacells_counts.csv, ss4_feature_metadata.csv, ss8_cell_metadata.csv, ss8_summarised_by_metacells_counts.csv, ss8_feature_metadata.csv, HH5_feature_metadata.csv, HH5_cell_metadata.csv, HH5_summarised_by_metacells_counts.csv, HH6_filtered_SEACells_integration_map.csv, HH6_ATAC_singlecell_integration_map.csv, HH6_seacells_seurat_integrated.RDS, ss8_filtered_SEACells_integration_map.csv, ss8_ATAC_singlecell_integration_map.csv, ss8_seacells_seurat_integrated.RDS, HH5_ATAC_singlecell_integration_map.csv, HH5_filtered_SEACells_integration_map.csv, HH5_seacells_seurat_integrated.RDS, HH7_ATAC_singlecell_integration_map.csv, HH7_filtered_SEACells_integration_map.csv, HH7_seacells_seurat_integrated.RDS, ss4_seacells_seurat_integrated.RDS, ss4_filtered_SEACells_integration_map.csv, ss4_ATAC_singlecell_integration_map.csv]]

    //combine all the summarised counts into one summarised counts file, check all feature metadata the same and write, combine all cell metadata too
    COMBINE_METACELL_COUNTS( ch_combined_input )
    
    // Filter peaks based on annotation and variability
    FILTER_PEAKS( COMBINE_METACELL_COUNTS.out )

    // Cluster peaks using Antler package
    CLUSTER_PEAKS( FILTER_PEAKS.out )

    CLUSTER_PEAKS.out.view()

    // Run homer motif enrichment on each peak module
    HOMER_MOTIF_ENRICHMENT( CLUSTER_PEAKS.out, ch_fasta )

    emit:
    clustered_peaks = CLUSTER_PEAKS.out
    
}
