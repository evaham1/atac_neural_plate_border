#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// calculate se on stages and full data
include {R as SE_CALCULATE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/calculate_se.R", checkIfExists: true) )

// compare intercluster variability between stages
include {R as COMPARE_INTERCLUSTER_VARIABILITY} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Measuring_variability/compare_intercluster_variability.R", checkIfExists: true) )

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// looking at differences in variability between stages

workflow COMPARE_VARIABILITY {
    take:
    input_ch

    main:
    // calculate se object for all stages + fulldata
    SE_CALCULATE( input_ch )

    // collect all se objects into a single channel
    ch_se_combined = SE_CALCULATE.out
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            .collect()
            .map { [[sample_id:'FullData'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    // compare variability/how many differential peaks we have at different stages
    COMPARE_INTERCLUSTER_VARIABILITY( ch_se_combined )
    
    emit:
    se_objects = SE_CALCULATE.out
}
