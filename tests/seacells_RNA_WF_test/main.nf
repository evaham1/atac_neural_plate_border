#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SEACELLS_RNA_WF } from '../../NF-downstream_analysis/subworkflows/local/PROCESSING/seacells_RNA_WF'
include { METADATA } from '../../NF-downstream_analysis/subworkflows/local/metadata'

Channel
    .value("../../NF-downstream_analysis/binary_knowledge_matrix_contam.csv")
    .set{ch_binary_knowledge_matrix}

workflow test_seacells_RNA_WF {

    METADATA( params.sample_sheet_test )
    ch_input = METADATA.out.metadata

    //ch_input.view()

    SEACELLS_RNA_WF( ch_input, ch_binary_knowledge_matrix )
}