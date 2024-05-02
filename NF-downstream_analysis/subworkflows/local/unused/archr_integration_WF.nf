#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// integration
//include {R as UNCON_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_unconstrained_integration.R", checkIfExists: true) )
// integration, cluster label and coaccessibility all in one script
include {R as INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_FULL_integration.R", checkIfExists: true) )

// checking label separation
//include {R as CLUSTER_IDENTIFY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )
include {R as DIM_RED_GENOMIC_SUBSETS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_dim_red_genomic_subsets.R", checkIfExists: true) )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow ARCHR_INTEGRATING_WF {
    take:
    input

    main:

    // Integrate full data and split stage data
    //UNCON_INTEGRATE ( input )
    INTEGRATE ( input )

    // Examine the resulting integration
    //CLUSTER_IDENTIFY ( UNCON_INTEGRATE.out ) // Visualise contributions of labels to each cluster and label clusters to summarise this

    // Try dimensionality reduction with different subsets of genome to see if they improve the separation of cell type labels
    //DIM_RED_GENOMIC_SUBSETS ( UNCON_INTEGRATE.out )

    //emit integrated ArchR objects:
    emit:
    integrated = INTEGRATE.out
}
