#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/1_preprocessing.R", checkIfExists: true) )
include {R as FILTERING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/2_filtering.R", checkIfExists: true) )
include {R as GENE_ACTIVITY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/3_gene_activity.R", checkIfExists: true) )

workflow PROCESSING {
    take:
    input

    main:
    input.view()

    PREPROCESSING( input )

    PREPROCESSING.out.view()

    FILTERING(PREPROCESSING.out)

    // filter input to only keep fragment files so can pass them to R processes that need them
    input
        .filter{ it[0].sample_id == 'NF-scATACseq_alignment_out' }
        .map {[it[0], it[1].collect{ file(it+"/outs/fragments.tsv.gz", checkIfExists: true) }]}
        .set {ch_fragments}

    ch_fragments.view()

    // combine the fragment files path with the rds objects from filtering
    // FILTERING.out
    //     .filter{ it[0].sample_id == 'NF-scATACseq_alignment_out' }
    //     .map {[it[0], it[1].collect{ file(it+"/outs/fragments.tsv.gz", checkIfExists: true) }]}
    //     .set {ch_filtering}

    GENE_ACTIVITY(ch_fragments)

    //emit:
    //preprocessing_out = PREPROCESSING.out
}
