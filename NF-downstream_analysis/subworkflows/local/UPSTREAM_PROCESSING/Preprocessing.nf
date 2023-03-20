#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PREPROCESS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_preprocessing.R", checkIfExists: true) )
include {EDIT_GTF} from "$baseDir/modules/local/edit_gtf/main"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    //emit full filtered and clustered dataset:
    emit:
    output = PREPROCESS.out
}