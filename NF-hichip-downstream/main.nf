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

// generate placeholder channel so can run generate_bins.R
Channel
    .value(params.gtf)
    .map { row -> [[sample_id:'dummy'], row] }
    .set{ch_dummy}


//
// WORKFLOW: Run main nf-core/NF-hichip-downstream analysis pipeline
//

workflow {

    //////////  Sample-generic bins generation and annotation  //////////

    // Turn gtf file into bed file - need to debug get bash syntax error!
    //GTF_TO_BED( ch_gtf )

    ch_dummy.view()

    // Generate bins
    GENERATE_BINS ( ch_dummy )

    GENERATE_BINS.out.view()
    // [[sample_id:WE_HiChip_r2], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/0e/f4d93e32d8bb18ca1c1dd1f781bc49/plots, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/0e/f4d93e32d8bb18ca1c1dd1f781bc49/rds_files]]
    // [[sample_id:WE_HiChip_r1], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/96/fc9346a0bd7dc603e92f2c45b79b56/plots, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/96/fc9346a0bd7dc603e92f2c45b79b56/rds_files]]
    // [[sample_id:WE_HiChip_r3], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/08/0bc2685c386e77728aa0c4ee7c90e3/plots, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/08/0bc2685c386e77728aa0c4ee7c90e3/rds_files]]
    // [[sample_id:NF_HiChip_r1], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/7e/72b5d662bcfce00a44041aa1181d27/plots, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/7e/72b5d662bcfce00a44041aa1181d27/rds_files]]
    // [[sample_id:NF_HiChip_r2], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/1c/edee56567009b8db00073ffbf28941/plots, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/1c/edee56567009b8db00073ffbf28941/rds_files]]
    // [[sample_id:NF_HiChip_r3], [/flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/8b/d96cb8c5093bde704e01eb344e5829/plots, /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/8b/d96cb8c5093bde704e01eb344e5829/rds_files]]

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
        // [[sample_id:NF_HiChip_r1], /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/cd/88cbc6d7ea5b20c6f52febaa5ec9fd/edited_ValidPairs.txt]
        // [[sample_id:WE_HiChip_r1], /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/c7/90d4da8bc07627e142e1dc2dabe34d/edited_ValidPairs.txt]
        // [[sample_id:WE_HiChip_r3], /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/10/9a2142476b94b5d081b47234204903/edited_ValidPairs.txt]
        // [[sample_id:WE_HiChip_r2], /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/cb/c9c6eb42cd96232528bc34b5de0c7d/edited_ValidPairs.txt]
        // [[sample_id:NF_HiChip_r2], /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/37/d9febe21517a7b599e57805d2ae542/edited_ValidPairs.txt]
        // [[sample_id:NF_HiChip_r3], /flask/scratch/briscoej/hamrude/atac_neural_plate_border/NF-hichip-downstream/work/22/321e6b00beda5a892912d42175ab62/edited_ValidPairs.txt]

    // Run HiCDCPlus on ValidPairs data to find significant interactions
    EDIT_VALIDPAIRS.out // EDIT_VALIDPAIRS.out: [[meta], [edited_ValidPairs.txt]]
        .combine( GENERATE_BINS.out )
        .map{[it[0], it[1] + it[2]]}
        .set { ch_validpairs }

    ch_validpairs.view()

    //here need some channel manipulation to combine GENERATE_BINS.out and EDIT_VALIDPAIRS.out



    //////////  Filter interesting interactions  //////////


}


/*
========================================================================================
    THE END
========================================================================================
*/
