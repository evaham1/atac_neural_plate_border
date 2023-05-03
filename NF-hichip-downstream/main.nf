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

// 2) Prep peak output and gtf and then intersect with bins to annotate them
include {GTF_TO_BED} from "$baseDir/modules/local/gtf_to_bed/main"
include {INTERSECT_BINS as INTERSECT_BINS_PEAKS} from "$baseDir/modules/local/intersect_bins/main"

// 3) Prep ValidPairs output from HiC pipeline and run with HiCDC+ to find significant interactions
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include {EDIT_VALIDPAIRS} from "$baseDir/modules/local/edit_ValidPairs/main"
include {R as LOOP_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/loop_calling_hicdc.R", checkIfExists: true) )

// 4) Filter bins to pick out interactions of interest

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

// generate placeholder channel so can run generate_bins.R
Channel
    .value(params.gtf)
    .map { row -> [[sample_id:'dummy'], row] }
    .set{ch_dummy}
//[[sample_id:dummy], /nemo/lab/briscoej/working/hamrude/raw_data/genomes/galgal6/tag_chroms.gtf]


//
// WORKFLOW: Run main nf-core/NF-hichip-downstream analysis pipeline
//

workflow {

    //////////  Sample-generic bins generation and annotation  //////////

    // Turn gtf file into bed file - need to debug get bash syntax error!
    //GTF_TO_BED( ch_gtf )

    // Generate bins
    GENERATE_BINS ( ch_dummy )
        // [[sample_id:dummy], [./work/c5/d2ca727dccb0dbb6645013d7e73c1e/plots, ./work/c5/d2ca727dccb0dbb6645013d7e73c1e/rds_files]]

    // Intersect bins with peaks
        //here need channel manipulation to combine peaks file from param and GENERATE_BINS.out
    GENERATE_BINS.out
        .map { row -> [row[0], row[1].findAll { it =~ ".*bed_files" }] }
        //.view() //[[sample_id:dummy], [bed_files]]
        .flatMap {it[1][0].listFiles()} //bins.bed
        .combine( ch_peaks )
        .view() //[bins.bed, FullData_PeakSet.bed]
        .set{ ch_peak_bins }

    INTERSECT_BINS_PEAKS( ch_peak_bins )

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

    // Run HiCDCPlus on ValidPairs data to find significant interactions
    EDIT_VALIDPAIRS.out // EDIT_VALIDPAIRS.out: [[sample_id:WE_HiChip_r1], edited_ValidPairs.txt]
        .map { row -> [row[0], [row[1]]] } //[[sample_id:WE_HiChip_r1], [edited_ValidPairs.txt]]
        .combine( GENERATE_BINS.out )
        .map{[it[0], it[1] + it[3]]}
        .set { ch_validpairs } //[[sample_id:NF_HiChip_r3], [edited_ValidPairs.txt, plots, rds_files]]

    //LOOP_CALL( ch_validpairs )

    //here need some channel manipulation to combine GENERATE_BINS.out and EDIT_VALIDPAIRS.out



    //////////  Filter interesting interactions  //////////


}


/*
========================================================================================
    THE END
========================================================================================
*/
