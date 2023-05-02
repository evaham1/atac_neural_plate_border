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
include {GTF_TO_BED} from "$baseDir/modules/local/gtf_to_bed/main"

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

workflow {

    //////////  Bins generation and annotation  //////////

    // Turn gtf file into bed file
    GTF_TO_BED( ch_gtf )

    // Generate bins
    GENERATE_BINS ( GTF_TO_BED.out )

    // Intersect bins with peaks
        //here need channel manipulation to combine peaks file from param and GENERATE_BINS.out

    // Intersect bins with genes
        //here need channel manipulation to combine GTF_TO_BED.out and GENERATE_BINS.out

    //////////  HiChip-sample specific analysis  //////////

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
        // [[sample_id:NF_HiChip_r1], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/89/191db4292187ae2b937f09cf34f809/edited_ValidPairs.txt, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/89/191db4292187ae2b937f09cf34f809/input]]
        // [[sample_id:WE_HiChip_r1], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/54/5cda609f53cfaa836a755b1f3f798d/edited_ValidPairs.txt, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/54/5cda609f53cfaa836a755b1f3f798d/input]]
        // [[sample_id:WE_HiChip_r2], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/3a/1edda7b589865c6ef67b033291421b/edited_ValidPairs.txt, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/3a/1edda7b589865c6ef67b033291421b/input]]
        // [[sample_id:NF_HiChip_r2], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/e6/84e3bec58c48e23a72869dbb67a110/edited_ValidPairs.txt, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/e6/84e3bec58c48e23a72869dbb67a110/input]]
        // [[sample_id:WE_HiChip_r3], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/c7/e9907b5bf54ae5559cb991235cdd0f/edited_ValidPairs.txt, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/c7/e9907b5bf54ae5559cb991235cdd0f/input]]
        // [[sample_id:NF_HiChip_r3], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/63/77a440049f7e59893b25d1f2b2f5d9/edited_ValidPairs.txt, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/63/77a440049f7e59893b25d1f2b2f5d9/input]]

    // Run HiCDCPlus on ValidPairs data to find significant interactions
        //here need some channel manipulation to combine GENERATE_BINS.out and EDIT_VALIDPAIRS.out



    //////////  Filter interesting interactions  //////////


}


/*
========================================================================================
    THE END
========================================================================================
*/
