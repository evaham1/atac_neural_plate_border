#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_split_stages.R", checkIfExists: true) )

include {R as CLUSTER_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )
include {R as PEAK_DIFF} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_diff_peaks.R", checkIfExists: true) )

// include {R as CLUSTER_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
// include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

// include {R as FILTER_CLUSTERS_1} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_clusters.R", checkIfExists: true) )
// include {R as FILTER_CLUSTERS_2} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_clusters.R", checkIfExists: true) )

// include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
// include {R as DOUBLETS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_doublets.R", checkIfExists: true) )


workflow QC_STAGES {
    take:
    input

    main:

    FILTER( input )

    ///     SPLIT STAGES    ///
    SPLIT_STAGES( FILTER.out )

    SPLIT_STAGES.out //[[meta], [plots, rds_files]]
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        //.view() //[[meta], [rds_files]]
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        .set { ch_split_stage }     
    
    ch_split_stage
       .view() //[[meta], Save-ArchR file]
    /////////////////////////////

    ///     CONFIRM IDENTITY OF LOW QUALITY CLUSTERS    ///
    CLUSTER_PREFILTER( ch_split_stage )
    GENE_SCORES( CLUSTER_PREFILTER.out )
    PEAK_CALL( CLUSTER_PREFILTER.out )
    PEAK_DIFF( PEAK_CALL.out )

    // /// ITERATIVELY CLUSTER AND FILTER STAGES ///
    // CLUSTER_PREFILTER( ch_split_stage )
    // FILTER_CLUSTERS_1( CLUSTER_PREFILTER.out ) // filtering round 1
    // CLUSTER_POSTFILTER( FILTER_CLUSTERS_1.out )
    // FILTER_CLUSTERS_2( CLUSTER_POSTFILTER.out ) // filtering round 2
    // CLUSTER( FILTER_CLUSTERS_2.out )
    // /////////////////////////


    



    // /// PLOT GENE SCORES ///
    // GENE_SCORES( CLUSTER.out )
    // /////////////////////////

    // /// IDENTIFY DOUBLETS ///
    // DOUBLETS( CLUSTER.out )
    // /////////////////////////

    // // extract rds objects
    // CLUSTER.out
    //     .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
    //     .flatMap {it[1][0].listFiles()}
    //     .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
    //     //.view() //CHECK THIS!
    //     .set {output_ch}

    //emit filtered and clustered stage objects:
    emit:
    //output = output_ch
    output = CLUSTER_PREFILTER.out
}
