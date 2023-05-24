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
include {SAMTOOLS_FAIDX} from "$baseDir/modules/local/samtools_faidx/main"
include {EXTRACT_PROMOTERS} from "$baseDir/modules/local/extract_promoters/main"
include {INTERSECT_BINS as INTERSECT_BINS_PEAKS} from "$baseDir/modules/local/intersect_bins/main"
include {INTERSECT_BINS as INTERSECT_BINS_PROMOTERS} from "$baseDir/modules/local/intersect_bins/main"

// 3) Prep ValidPairs output from HiC pipeline and run with HiCDC+ to find differential significant interactions
include { METADATA } from "$baseDir/subworkflows/local/metadata"
include {EDIT_VALIDPAIRS} from "$baseDir/modules/local/edit_ValidPairs/main"
include {R as LOOP_CALL} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/loop_calling_hicdc.R", checkIfExists: true) )
include {R as DIFF_LOOPS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/diff_loops_hicdc.R", checkIfExists: true) )

// 4) Filter bins to pick out interactions of interest
include {R as INVESTIGATE_LOOPS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/investigate_loops.R", checkIfExists: true) )


//
// SET CHANNELS
//

// set channel for genome index file
Channel
    .value(params.index)
    .set{ch_index}

// set channel for genome gtf file
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

    // Extract promoters of genes from gtf and turn into bed file
    ch_gtf
        .combine( ch_index )
        //.view() // [/nemo/lab/briscoej/working/hamrude/raw_data/genomes/galgal6/tag_chroms.gtf, /nemo/lab/briscoej/working/hamrude/raw_data/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa.fai]
        .set{ ch_extract_bins }
    EXTRACT_PROMOTERS( ch_extract_bins )

    // Generate bins
    GENERATE_BINS ( ch_dummy )
        // [[sample_id:dummy], [bed_files, plots, rds_files]]

    // Extract bins bed file from generate_bins.R output
    GENERATE_BINS.out
        .map { row -> [row[0], row[1].findAll { it =~ ".*bed_files" }] }
        //.view() //[[sample_id:dummy], [bed_files]]
        .flatMap {it[1][0].listFiles()} //bins.bed
        .set{ ch_bins }

    // Intersect bins with peaks
    ch_bins
        .combine( ch_peaks )
        //.view() //[bins.bed, FullData_PeakSet.bed]
        .set{ ch_peak_bins }
    INTERSECT_BINS_PEAKS( ch_peak_bins )

    // Intersect bins with promoters
    ch_bins
        .combine( EXTRACT_PROMOTERS.out )
        //.view() //[bins.bed, tag_chroms.bed]
        .set{ ch_promoters_bins }
    INTERSECT_BINS_PROMOTERS( ch_promoters_bins )

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
        // [[sample_id:WE_HiChip_r2], f6/5f387e28158d20d54bdb9825be3137/WE_HiChip_r2_edited.allValidPairs]
        // [[sample_id:WE_HiChip_r1], 63/5c215c3285e67ba1b10da286aaf5a2/WE_HiChip_r1_edited.allValidPairs]
        // [[sample_id:WE_HiChip_r3], b6/507a4d889f71a341d05dadbc0c5de8/WE_HiChip_r3_edited.allValidPairs]
        // [[sample_id:NF_HiChip_r1], ca/387797417c462833e88dbc97965220/NF_HiChip_r1_edited.allValidPairs]
        // [[sample_id:NF_HiChip_r2], 19/4723cf3020c3fda9de571a41609486/NF_HiChip_r2_edited.allValidPairs]

    // Run HiCDCPlus on ValidPairs data to find significant interactions
    EDIT_VALIDPAIRS.out // EDIT_VALIDPAIRS.out: [[sample_id:WE_HiChip_r2], WE_HiChip_r2_edited.allValidPairs]
        .map { row -> [row[0], [row[1]]] } //[[sample_id:WE_HiChip_r2], [WE_HiChip_r2_edited.allValidPairs]]
        .combine( GENERATE_BINS.out ) //[[sample_id:dummy], [bed_files, plots, rds_files]]
        .map{[it[0], it[1] + it[3]]}
        .set { ch_validpairs } //[[sample_id:WE_HiChip_r1], [WE_HiChip_r2_edited.allValidPairs, bed_files, plots, rds_files]]

    LOOP_CALL( ch_validpairs )
    LOOP_CALL.out.view()

    // //////////  Find differential interactions between NF and WE  //////////

    LOOP_CALL.out
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0][1]}
        .collect()
        .map { [[sample_id:'AllSamples'], it] } //
        .view()
        .set{ ch_interactions_combined }

    //[[sample_id:AllSamples], [WE_HiChip_r1_HiCDC_output_filtered.txt, NF_HiChip_r1_HiCDC_output.txt.gz, WE_HiChip_r3_HiCDC_output.txt.gz, NF_HiChip_r2_HiCDC_output_filtered.txt, NF_HiChip_r3_HiCDC_output.txt.gz, WE_HiChip_r2_HiCDC_output_filtered.txt]]
    DIFF_LOOPS( ch_interactions_combined )

    //////////  Pull out interesting interactions  //////////

    // DIFF_LOOPS.out
    //     .map { row -> [row[0], row[1].findAll { it =~ ".*rds_files" }[0]] } //[[sample_id:WE_HiChip_r2], rds_files]
    //     .combine( INTERSECT_BINS_PEAKS.out )
    //     .combine( INTERSECT_BINS_PROMOTERS.out )
    //     .combine( ch_bins )
    //     .map{ [ it[0], [it[1], it[2], it[3], it[4]] ] }
    //     .view() // [[sample_id:WE_HiChip_r1], [rds_files, FullData_PeakSet_bins_intersected.bed, tag_chroms_bins_intersected.bed, bins.bed]]
    //     .set{ ch_intersect }

    // INVESTIGATE_LOOPS( ch_intersect )

}

/*
========================================================================================
    THE END
========================================================================================
*/
