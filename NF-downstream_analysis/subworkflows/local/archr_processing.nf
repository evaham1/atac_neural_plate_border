#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as ARCHR_PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_ArchR_preprocessing.R", checkIfExists: true) )
include {R as ARCHR_FILTERING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/2_ArchR_filtering.R", checkIfExists: true) )
include {R as ARCHR_VISUALISE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/3_ArchR_visualisations.R", checkIfExists: true) )
include {EDIT_GTF} from "$baseDir/modules/local/edit_gtf/main"

workflow ARCHR_PROCESSING {
    take:
    input

    main:
    input // [[meta], [cellranger_output, galgal_reference]]
        .set {ch_input}

    EDIT_GTF ( input ) //edits the gtf file to add 'chr' to chromosome names

    EDIT_GTF.out 
        .combine(ch_input) //[[meta], temp.gtf, [meta], cellranger_output, galgal_reference]]
        .map{[it[0], it[[1]] + it[3]]} //[[meta], [temp.gtf, cellranger_output, galgal_reference]]
        .set {ch_input_modified} // ch_metadata: [[meta], [cellranger_output, gtf]]

    ARCHR_PREPROCESSING( ch_input_modified )
    ARCHR_FILTERING( ARCHR_PREPROCESSING.out )
    ARCHR_VISUALISE( ARCHR_FILTERING.out )

    //emit:
    //signac_predicted_gex = GEX_FILTERING.out
}
