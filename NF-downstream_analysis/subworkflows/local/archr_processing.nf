#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {EDIT_GTF} from "$baseDir/modules/local/edit_gtf/main"

include {R as ARCHR_PREPROCESSING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_preprocessing.R", checkIfExists: true) )
include {R as ARCHR_DOUBLETS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_doublets.R", checkIfExists: true) )
include {R as ARCHR_FILTERING} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filtering.R", checkIfExists: true) )
include {R as ARCHR_DOUBLETS_FILTERED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_doublets.R", checkIfExists: true) )

include {R as ARCHR_CLUSTERING_PREFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as ARCHR_CLUSTERING_POSTFILTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
include {R as ARCHR_CLUSTERING_POSTFILTER_TWICE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

include {R as ARCHR_FILTER_CLUSTERS_1} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_clusters.R", checkIfExists: true) )
include {R as ARCHR_FILTER_CLUSTERS_2} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_filter_clusters.R", checkIfExists: true) )

include {R as ARCHR_GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )


workflow WF_ARCHR_PROCESSING {
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

    // creates arrow files and ArchR project filtered with generous thresholds
    ARCHR_PREPROCESSING( ch_input_modified )
    // plots whole sample QC metrics (+ add filtering?)
    ARCHR_FILTERING( ARCHR_PREPROCESSING.out )

    // iterative clustering and filtering poor quality clusters
    ARCHR_CLUSTERING_PREFILTER( ARCHR_FILTERING.out )
    ARCHR_FILTER_CLUSTERS_1( ARCHR_CLUSTERING_PREFILTER.out ) // filtering round 1
    ARCHR_CLUSTERING_POSTFILTER( ARCHR_FILTER_CLUSTERS_1.out )
    ARCHR_FILTER_CLUSTERS_2( ARCHR_CLUSTERING_POSTFILTER.out ) // filtering round 2
    ARCHR_CLUSTERING_POSTFILTER_TWICE( ARCHR_FILTER_CLUSTERS_2.out )
    ARCHR_DOUBLETS_FILTERED( ARCHR_FILTER_CLUSTERS_2.out ) // see if adding doublet scores after filtering any better

    // plots using gene scores
    ARCHR_GENE_SCORES( ARCHR_CLUSTERING_POSTFILTER.out )

    // extract rds objects
      ARCHR_CLUSTERING_POSTFILTER.out //[[sample_id:NF-scATACseq_alignment_out], [../ArchRLogs, ../Rplots.pdf, ../rds_files]]
        .map {row -> [row[0], row[1].findAll { it =~ ".*rds_files" }]} //[[sample_id:NF-scATACseq_alignment_out], [../rds_files]]
        .flatMap {it[1][0].listFiles()}
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        //.view() //[[sample_id:FullData], /rds_files/FullData_Save-ArchR]
        .set {output_ch}

    //emit full filtered and clustered dataset:
    emit:
    output = output_ch
}