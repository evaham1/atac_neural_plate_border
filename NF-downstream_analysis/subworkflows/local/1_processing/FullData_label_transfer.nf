#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//include {R as TRANSFER_LABELS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )
//include {R as SUBSET} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

include {R as CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )

include {R as GENE_SCORES} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_gene_scores.R", checkIfExists: true) )
include {R as DOUBLETS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_doublets.R", checkIfExists: true) )


workflow FULL_PROCESSING {
    take:
    input

    main:

    /// TRANSFER LABELS FROM STAGES TO FILTER FULL DATA ///
    // transfer labels

    // subset script
    /////////////////////////

    /// CLUSTER ///
    CLUSTER( FILTER_CLUSTERS_2.out )
    /////////////////////////
    
    /// PLOT GENE SCORES ///
    GENE_SCORES( CLUSTER.out )
    /////////////////////////

    /// IDENTIFY DOUBLETS ///
    DOUBLETS( CLUSTER.out )
    /////////////////////////

    // extract rds objects
    CLUSTER.out //[[sample_id:NF-scATACseq_alignment_out], [../ArchRLogs, ../Rplots.pdf, ../rds_files]]
        .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }] } //[[sample_id:NF-scATACseq_alignment_out], [../rds_files]]
        .flatMap { it[1][0].listFiles() }
        .map { row -> [[sample_id:row.name.replaceFirst(~/_[^_]+$/, '')], row] }
        //.view() //[[sample_id:FullData], /rds_files/FullData_Save-ArchR]
        .set { output_ch }

    //emit filtered and clustered stage objects:
    emit:
    output = output_ch
}
