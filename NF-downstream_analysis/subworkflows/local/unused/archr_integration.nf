#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// integration
include {R as UNCON_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_unconstrained_integration.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )

// remove contamination
include {R as SUBSET_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_subsetting.R", checkIfExists: true) )

// rerun clustering and check integration on subset
include {R as CLUSTER_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY_FILTERED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow INTEGRATING {
    take:
    input

    main:

    // Integrate full data and split stage data
    UNCON_INTEGRATE ( input )

    // Examine the resulting integration
    CLUSTER_IDENTIFY ( UNCON_INTEGRATE.out ) // Visualise contributions of labels to each cluster and label clusters to summarise this
    
////////////    FILTER OUT CONTAMINATION    ///////////////////////

    // Filter contaminating cells from all channels and re-cluster all channels
    SUBSET_INTEGRATION ( UNCON_INTEGRATE.out )
    CLUSTER_INTEGRATION ( SUBSET_INTEGRATION.out )

    // Examine the resulting integration
    CLUSTER_IDENTIFY_FILTERED ( CLUSTER_INTEGRATION.out ) // Visualise contributions of labels to each cluster and label clusters to summarise this

    //emit integrated ArchR objects:
    emit:
    integrated = UNCON_INTEGRATE.out
    integrated_filtered = CLUSTER_IDENTIFY_FILTERED.out
}