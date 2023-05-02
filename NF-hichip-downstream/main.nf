#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/NF-hichip-downstream
========================================================================================
    Github : 
    Website: 
    Slack  : 
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

// 1) Create bins
include {R as GENERATE_BINS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/generate_bins.R", checkIfExists: true) )

// 2) Prep ValidPairs output from HiC pipeline and run with HiCDC+ to find significant interactions
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include {EDIT_VALIDPAIRS} from "$baseDir/modules/local/edit_ValidPairs/main"

// 2) Prep peak output and gtf and then intersect with bins to annotate them

// 3) Filter bins to pick out interactions of interest

//
// SET CHANNELS
//

// set channel for gtf file
Channel
    .value(params.gtf)
    .set{ch_gtf}

// set channel for peaks.bed file
Channel
    .value(params.peaks)
    .set{ch_peaks}


//
// WORKFLOW: Run main nf-core/NF-hichip-downstream analysis pipeline
//

workflow A {

    // Read in ValidPairs data
    METADATA( params.sample_sheet_validpairs )
        // [[sample_id:WE_HiChip_r1], [/flask/scratch/briscoej/thierya/atac_neural_plate_border/output/NF-hichip_alignment/hicpro/valid_pairs/WE_HiChip_r1.allValidPairs]]
        // [[sample_id:WE_HiChip_r2], [/flask/scratch/briscoej/thierya/atac_neural_plate_border/output/NF-hichip_alignment/hicpro/valid_pairs/WE_HiChip_r2.allValidPairs]]
        // [[sample_id:WE_HiChip_r3], [/flask/scratch/briscoej/thierya/atac_neural_plate_border/output/NF-hichip_alignment/hicpro/valid_pairs/WE_HiChip_r3.allValidPairs]]
        // [[sample_id:NF_HiChip_r1], [/flask/scratch/briscoej/thierya/atac_neural_plate_border/output/NF-hichip_alignment/hicpro/valid_pairs/NF_HiChip_r1.allValidPairs]]
        // [[sample_id:NF_HiChip_r2], [/flask/scratch/briscoej/thierya/atac_neural_plate_border/output/NF-hichip_alignment/hicpro/valid_pairs/NF_HiChip_r2.allValidPairs]]
        // [[sample_id:NF_HiChip_r3], [/flask/scratch/briscoej/thierya/atac_neural_plate_border/output/NF-hichip_alignment/hicpro/valid_pairs/NF_HiChip_r3.allValidPairs]]

    // Edit ValidPairs data to add 'chr' to chromosome names
    EDIT_VALIDPAIRS ( METADATA.out )

    // Generate bins
    GENERATE_BINS ()

}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    A ()
}


/*
========================================================================================
    THE END
========================================================================================
*/
