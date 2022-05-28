#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_split_stages.R", checkIfExists: true) )

include {R as FILTER_CLUSTER_LOOP} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_cluster_filter_loop.R", checkIfExists: true) )

include {R as CLUSTER_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as PEAK_CALL_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as HEATMAP_PEAKS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/plot_marker_heatmaps.R", checkIfExists: true) )
include {R as HEATMAP_GEX} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/plot_marker_heatmaps.R", checkIfExists: true) )


workflow FILTERING {
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

    ///     FILTER CLUSTERS     ///
    FILTER_CLUSTER_LOOP( ch_split_stage )

    ///     CHECK QUALITY     ///
    CLUSTER_POSTFILTER( FILTER_CLUSTER_LOOP.out )
    GENE_SCORES_POSTFILTER( CLUSTER_POSTFILTER.out )
    HEATMAP_GEX( PEAK_CALL_POSTFILTER.out )

    PEAK_CALL_POSTFILTER( CLUSTER_POSTFILTER.out )
    HEATMAP_PEAKS( PEAK_CALL_POSTFILTER.out )

    emit:
    output = PEAK_CALL_POSTFILTER.out

}