#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as SPLIT_ATAC} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/split_seurat.R", checkIfExists: true) )
include {R as SPLIT_RNA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/split_seurat.R", checkIfExists: true) )
include {R as INTEGRATE_RNA} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/4_integrate_rna.R", checkIfExists: true) )

workflow INTEGRATE_SPLIT_PROCESS {
    take:
    ch_atac              //Channel: [seurat_atac]
    ch_rna              //Channel: [seurat_rna]
    ch_cellranger        //Channel: [cellranger_output]

    main:

    // adding metadata
    ch_atac
        .map{row -> [[sample_id:'atac'], row]} 
        .set{ ch_atac_input } // [[sample_id:atac], seurat_GeneActivity.RDS]
    ch_rna
        .map{row -> [[sample_id:'rna'], row]} 
        .set{ ch_rna_input } // [[sample_id:rna], seurat_label_transfer.RDS]

    // split the atac and rna objects
    SPLIT_ATAC ( ch_atac_input )
    SPLIT_RNA ( ch_rna_input )

    // change metadata so its the stage

    SPLIT_ATAC.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_atac_out } 
    SPLIT_RNA.out
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]}
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/\.[^\.]+$/, '')], row] }
        .set { ch_rna_out } 

    // combine split rna and atac objects together and with cellranger output
    ch_atac_out
        .combine(ch_rna_out, by:0)
        .combine(ch_cellranger)
        .map{[it[0], it[1] + it[2] + it[3] ]}
        .set { ch_integrate_input } 

    ch_integrate_input.view()
    // run integration on each stage individually
    INTEGRATE_RNA ( ch_integrate_input )

    emit:
    signac_integrated = INTEGRATE_RNA.out
}