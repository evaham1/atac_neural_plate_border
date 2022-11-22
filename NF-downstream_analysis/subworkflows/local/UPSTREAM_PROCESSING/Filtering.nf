#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as FILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_split_stages.R", checkIfExists: true) )

include {R as FILTER_CLUSTER_LOOP} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_cluster_filter_loop.R", checkIfExists: true) )
include {R as FILTER_FULL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_full_data.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow FILTERING {
    take:
    input

    main:

    input
        .set {ch_input}

    ///     FILTER nFRAGS    ///
    FILTER( input )

    ///     SPLIT STAGES    ///
    SPLIT_STAGES( FILTER.out )
    SPLIT_STAGES.out //[[meta], [plots, rds_files]]
        .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] }
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        .set { ch_split_stage }     

    ///     FILTER CLUSTERS IN INDIVIDUAL STAGES     ///
    FILTER_CLUSTER_LOOP( ch_split_stage )

    ///     FILTER FULL DATA USING CELL IDS FROM STAGES     ///
    ch_combined = FILTER_CLUSTER_LOOP.out // Collect rds files from all stages
            .concat(ch_input)
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            .collect()
            .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    FILTER_FULL ( ch_combined ) // filter full data

    ch_output = FILTER_CLUSTER_LOOP.out // Collect rds files from all stages
        .concat(FILTER_FULL.out)
        .map{[it[0], it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]]} //[ [[meta: HH5], ATAC.rds] , [[meta: HH6], ATAC.rds], [[meta: FullData], ATAC.rds]]
        .map{ [ it[0], it[[1]] ] }
        .view()

    emit:
    output = ch_output

}
