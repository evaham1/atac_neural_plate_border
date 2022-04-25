#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {R as SPLIT_STAGES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_split_stages.R", checkIfExists: true) )
include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )


workflow STAGE_PROCESSING {
    take:
    input

    main:
    SPLIT_STAGES( input )

    SPLIT_STAGES.out //[[meta], [plots, rds_files]]
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        //.view() //[[meta], [rds_files]]
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        .set { ch_split_stage }     
    
    ch_split_stage
       .view() //[[meta], Save-ArchR file]

    // cluster individual stages
    CLUSTER( ch_split_stage )
    
    // gene score plots for individual stages
    GENE_SCORES( CLUSTER.out )

    // extract rds objects
    CLUSTER.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        //.view() //CHECK THIS!
        .set {output_ch}

    //emit full filtered and clustered dataset:
    emit:
    output = output_ch
}
