#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// R scripts to run differential expression tests to find NC and PPR-specific enhancers
include {R as CALCULATE_SE} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Peak_calling/calculate_se.R", checkIfExists: true) )
include {R as FIND_ENHANCERS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Finding_enhancers/finding_enhancers.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow FIND_ENHANCERS_WF {
    take:
    input

    main:
    
    CALCULATE_SE( input )
    FIND_ENHANCERS( CALCULATE_SE.out )

    emit:
    se_object = CALCULATE_SE.out
    enhancers_output = FIND_ENHANCERS.out

}
