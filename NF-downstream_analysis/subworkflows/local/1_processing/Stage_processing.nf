#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_split_stages.R", checkIfExists: true) )

include {R as CLUSTER_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as PEAK_CALL_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as PEAK_DIFF_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_diff_peaks.R", checkIfExists: true) )

include {R as FILTER_CLUSTERS_1} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_clusters.R", checkIfExists: true) )
include {R as CLUSTER_PREFILTER_1} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as FILTER_CLUSTERS_2} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_clusters.R", checkIfExists: true) )

include {R as CLUSTER_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as PEAK_CALL_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as PEAK_DIFF_POSTFILTER} from "$baseDir/modules/local/r/main"                addParams(script: file("$baseDir/bin/Peak_calling/ArchR_diff_peaks.R", checkIfExists: true) )

workflow QC_STAGES {
    take:
    input

    main:

    ///     FILTER nFRAGS    ///
    FILTER( input )

    ///     SPLIT STAGES    ///
    SPLIT_STAGES( FILTER.out )
    SPLIT_STAGES.out //[[meta], [plots, rds_files]]
        .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        .set { ch_split_stage }     

    ///     CONFIRM IDENTITY OF LOW QUALITY CLUSTERS    ///
    CLUSTER_PREFILTER( ch_split_stage )
    //GENE_SCORES_PREFILTER( CLUSTER_PREFILTER.out )
    //PEAK_CALL_PREFILTER( CLUSTER_PREFILTER.out )
    //PEAK_DIFF_PREFILTER( PEAK_CALL_PREFILTER.out )

    ///     FILTER CLUSTERS TWICE     ///
    FILTER_CLUSTERS_1( CLUSTER_PREFILTER.out )
    CLUSTER_PREFILTER_1( FILTER_CLUSTERS_1.out )
    FILTER_CLUSTERS_2( CLUSTER_PREFILTER_1.out )

    ///     PLOTS FOR FILTERED DATA    ///
    CLUSTER_POSTFILTER( FILTER_CLUSTERS_2.out )
    //GENE_SCORES_POSTFILTER( CLUSTER_POSTFILTER.out )
    //PEAK_CALL_POSTFILTER( CLUSTER_POSTFILTER.out )
    //PEAK_DIFF_POSTFILTER( PEAK_CALL_POSTFILTER.out )

    // emit filtered and clustered stage objects:
    emit:
    //unfiltered_output = PEAK_CALL_PREFILTER.out
    //unfiltered_gex = GENE_SCORES_PREFILTER.out
    //unfiltered_differential_peaks = PEAK_DIFF_PREFILTER.out
    //output = PEAK_CALL_POSTFILTER.out
    //gex = GENE_SCORES_POSTFILTER.out
    //differential_peaks = PEAK_DIFF_POSTFILTER.out
    output = CLUSTER_POSTFILTER.out
}
