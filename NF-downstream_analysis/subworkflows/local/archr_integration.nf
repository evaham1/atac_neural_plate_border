#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// integration
include {R as UNCON_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_unconstrained_integration.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )
include {R as INTEGRATION_CLUSTERS_COMPARE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_compare_clusters_and_labels.R", checkIfExists: true) )

// remove contamination
include {R as SUBSET_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_subsetting.R", checkIfExists: true) )
include {R as CLUSTER_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_utilities/ArchR_clustering.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY_FILTERED} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow INTEGRATING {
    take:
    input

    main:

    // Integrate full data and split stage data
    UNCON_INTEGRATE ( input )

    // Examine the resulting integration
    CLUSTER_IDENTIFY ( UNCON_INTEGRATE.out ) // Visualise contributions of labels to each cluster and label clusters to summarise this
    INTEGRATION_CLUSTERS_COMPARE ( UNCON_INTEGRATE.out ) // Examine the relationship between clusters and labels
    
////////////    FILTER OUT CONTAMINATION    ///////////////////////

    // Filter contaminating cells from all channels and re-cluster all channels
    SUBSET_INTEGRATION ( UNCON_INTEGRATE.out )
    CLUSTER_INTEGRATION ( SUBSET_INTEGRATION.out )

    // Examine the resulting integration
    CLUSTER_IDENTIFY ( CLUSTER_INTEGRATION.out ) // Visualise contributions of labels to each cluster and label clusters to summarise this
    INTEGRATION_CLUSTERS_COMPARE ( CLUSTER_INTEGRATION.out ) // Examine the relationship between clusters and labels
    
    //emit integrated ArchR objects:
    emit:
    integrated = UNCON_INTEGRATE.out
    integrated_filtered = CLUSTER_INTEGRATION.out
}
