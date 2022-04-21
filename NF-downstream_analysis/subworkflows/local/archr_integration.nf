#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as UNCON_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_unconstrained_integration.R", checkIfExists: true) )
include {R as CLUSTER_IDENTIFY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_cluster_identities.R", checkIfExists: true) )

include {R as SUBSET_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
include {R as CLUSTER_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )


workflow INTEGRATING {
    take:
    input_ch

    main:
    // Integrate full data and split stage data
    UNCON_INTEGRATE ( input_ch )

    // Label clusters based on most frequent label within each cluster
    CLUSTER_IDENTIFY ( UNCON_INTEGRATE.out )
    
    // Filter contaminating cells from all channels and re-cluster all channels
    SUBSET_INTEGRATION ( INTEGRATE.out )
    CLUSTER_INTEGRATION ( SUBSET_INTEGRATION.out )

    //emit integrated ArchR objects:
    emit:
    archr_integrated_full = CLUSTER_INTEGRATION.out
}
