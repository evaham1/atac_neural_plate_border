#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as ARCHR_INTEGRATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/ArchR_integration.R", checkIfExists: true) )

workflow ARCHR_INTEGRATION {
    take:
    ATAC_list_data //Channel: [[meta], [HH5, rds_dir]]
    RNA_list_data //Channel: [[meta], [plot_dir, rds_dir]]

    main:
    // need to add here adding full datasets to both ATAC and RNA lists

    // combine ATAC list with RNA list one by one
    ATAC_list_data
        .view()

    ARCHR_INTEGRATE ( ATAC_list_data )
    
    // script to integrate datasets, run on full data and on subsets


    //emit integrated ArchR objects:
    emit:
    archr_integrated_full = ARCHR_INTEGRATE.out
}
