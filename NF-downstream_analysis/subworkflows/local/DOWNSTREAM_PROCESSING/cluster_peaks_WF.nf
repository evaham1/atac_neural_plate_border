#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R script to combine metacell outputs
include {R as COMBINE_METACELL_COUNTS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Metacell_processes/Combine_summarised_counts.R", checkIfExists: true) )

// R scripts to filter peaks and then cluster them using Antler
include {R as FILTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/cluster_peaks/1_filter_peaks.R", checkIfExists: true) )
include {R as CLUSTER_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/cluster_peaks/2_antler_calculate_peak_modules.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// cluster peaks

workflow CLUSTER_PEAKS_WF {
    take:
    input //should be SEACELLS_ATAC_WF.out.seacell_outputs_named

    main:

    //input.view()
        // [[sample_id:HH7], 58/b3df47b5a798e1acaa3666df878bdb/csv_files]
        // [[sample_id:HH6], 41/8c7ec3c7a7c3bb79dcd52e4bad02b9/csv_files]
        // [[sample_id:ss8], a8/f7a307efa1093759a27fd61ea09350/csv_files]
        // [[sample_id:ss4], 71/b91188feb85928e3c2096811627513/csv_files]
        // [[sample_id:HH5], b5/52af06cc3559a5d45a976102bc509f/csv_files]

    // take the exported_data outputs from SEACell_computation of the ATAC
    ch_metacells_combined = input // Collect csv files from all stages
        .map{ meta, data -> data.listFiles() } //list files inside only folder in output (csv_files)
            // [HH6_feature_metadata.csv, HH6_summarised_by_metacells_counts.csv, HH6_cell_metadata.csv]
            // [ss4_cell_metadata.csv, ss4_summarised_by_metacells_counts.csv, ss4_feature_metadata.csv]
            // [ss8_cell_metadata.csv, ss8_summarised_by_metacells_counts.csv, ss8_feature_metadata.csv]
            // [HH7_feature_metadata.csv, HH7_summarised_by_metacells_counts.csv, HH7_cell_metadata.csv]
            // [HH5_feature_metadata.csv, HH5_cell_metadata.csv, HH5_summarised_by_metacells_counts.csv]
        .flatten() //removes square brackets from each array
        .collect() // puts all arrays together
        .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    //ch_metacells_combined.view()
        //[[sample_id:FullData], [HH7_feature_metadata.csv, HH7_summarised_by_metacells_counts.csv, HH7_cell_metadata.csv, HH6_feature_metadata.csv, HH6_summarised_by_metacells_counts.csv, HH6_cell_metadata.csv, HH5_feature_metadata.csv, HH5_cell_metadata.csv, HH5_summarised_by_metacells_counts.csv, ss8_cell_metadata.csv, ss8_summarised_by_metacells_counts.csv, ss8_feature_metadata.csv, ss4_cell_metadata.csv, ss4_summarised_by_metacells_counts.csv, ss4_feature_metadata.csv]]

    //combine all the summarised counts into one summarised counts file, check all feature metadata the same and write, combine all cell metadata too
    COMBINE_METACELL_COUNTS( ch_metacells_combined )
    
    // Filter peaks based on annotation and variability
    FILTER_PEAKS( input )

    // Cluster peaks using Antler package
    CLUSTER_PEAKS( FILTER_PEAKS.out )

    emit:
    test_output = input
    clustered_peaks = CLUSTER_PEAKS.out
    
}
