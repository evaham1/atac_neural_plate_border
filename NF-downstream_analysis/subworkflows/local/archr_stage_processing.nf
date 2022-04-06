#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {R as ARCHR_SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_split_stages.R", checkIfExists: true) )
include {R as ARCHR_CLUSTERING_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as ARCHR_GENE_SCORES_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )


workflow ARCHR_STAGE_PROCESSING {
    take:
    input

    main:
    ARCHR_SPLIT_STAGES( input )

    ARCHR_SPLIT_STAGES.out //[[meta], [plots, rds_files]]
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        //.view() //[[meta], [rds_files]]
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        .set { ch_split_stage }     
    
    //ch_split_stage
    //    .view() //[[meta], Save-ArchR file]

    // cluster individual stages
    ARCHR_CLUSTERING_STAGES( ch_split_stage )
    
    // gene score plots for individual stages
    ARCHR_GENE_SCORES_STAGES( ARCHR_CLUSTERING_STAGES.out )

    // extract rds objects
    ARCHR_CLUSTERING_STAGES.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        .view() //CHECK THIS!
        .set {output_ch}

    //emit full filtered and clustered dataset:
    emit:
    output = output_ch
}
