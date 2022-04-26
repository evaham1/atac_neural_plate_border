#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {EDIT_GTF} from "$baseDir/modules/local/edit_gtf/main"

include {R as PREPROCESS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_preprocessing.R", checkIfExists: true) )
include {R as FILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )

include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )


workflow PREPROCESSING {
    take:
    input

    main:

    /// SET UP INPUT DATA ///
    input // [[meta], [cellranger_output, galgal_reference]]
        .set {ch_input}

    EDIT_GTF ( input ) //edits the gtf file to add 'chr' to chromosome names

    EDIT_GTF.out 
        .combine(ch_input) //[[meta], temp.gtf, [meta], cellranger_output, galgal_reference]]
        .map{[it[0], it[[1]] + it[3]]} //[[meta], [temp.gtf, cellranger_output, galgal_reference]]
        .set {ch_input_modified} // ch_metadata: [[meta], [cellranger_output, gtf]]
    /////////////////////////

    /// CREATE ARCHR PROJECT ///
    PREPROCESS( ch_input_modified )
    /////////////////////////

    /// FILTER ON ALL SAMPLES ///
    FILTER( PREPROCESS.out )
    /////////////////////////
    
    // DOUBLETS_FILTERED( FILTER_CLUSTERS_2.out ) // see if adding doublet scores after filtering any better

    // // plots using gene scores
    // GENE_SCORES( CLUSTER.out )

    // // extract rds objects
    //   CLUSTER.out //[[sample_id:NF-scATACseq_alignment_out], [../ArchRLogs, ../Rplots.pdf, ../rds_files]]
    //     .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] } //[[sample_id:NF-scATACseq_alignment_out], [../rds_files]]
    //     .flatMap { it[1][0].listFiles() }
    //     .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
    //     //.view() //[[sample_id:FullData], /rds_files/FullData_Save-ArchR]
    //     .set { output_ch }

    //emit full filtered and clustered dataset:
    emit:
    output = FILTER.out
}