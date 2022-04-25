#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {R as PEAK_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/ArchR_peak_calling.R", checkIfExists: true) )

workflow PEAK_CALLING {
    take:
    input_ch

    main:
    // run peak calling on clusters
    PEAK_CALL ( input_ch )
    
    
    emit:
    archr_peaks = PEAK_CALL.out
}
