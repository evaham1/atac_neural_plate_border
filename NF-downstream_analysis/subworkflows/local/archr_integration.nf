#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_integration.R", checkIfExists: true) )
include {R as SUBSET_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
include {R as CLUSTER_INTEGRATION} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )


workflow INTEGRATING {
    take:
    input_ch

    main:
    // Integrate full data and split stage data
    INTEGRATE ( input_ch )
    
    // Filter contaminating cells from all channels and re-cluster all channels
    SUBSET_INTEGRATION ( INTEGRATE.out )
    CLUSTER_INTEGRATION ( SUBSET_INTEGRATION.out )

    //emit integrated ArchR objects:
    emit:
    archr_integrated_full = CLUSTER_INTEGRATION.out
}
