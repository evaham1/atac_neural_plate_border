#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as ARCHR_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_integration.R", checkIfExists: true) )
include {R as ARCHR_INTEGRATE_SUBSET} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_subsetting.R", checkIfExists: true) )
include {R as ARCHR_INTEGRATE_CLUSTER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/ArchR_preprocessing/ArchR_clustering.R", checkIfExists: true) )


workflow ARCHR_INTEGRATION {
    take:
    input_ch

    main:
    // Integrate full data and split stage data
    ARCHR_INTEGRATE ( input_ch )

    ARCHR_INTEGRATE.out.view()
    
    // Filter contaminating cells from all channels and re-cluster all channels
    ARCHR_INTEGRATE_SUBSET ( ARCHR_INTEGRATE.out )
    ARCHR_INTEGRATE_CLUSTER ( ARCHR_INTEGRATE_SUBSET.out )

    //emit integrated ArchR objects:
    emit:
    archr_integrated_full = ARCHR_INTEGRATE_CLUSTER.out
}
