#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/downstream
========================================================================================
    Github : https://github.com/nf-core/downstream
    Website: https://nf-co.re/downstream
    Slack  : https://nfcore.slack.com/channels/downstream
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { METADATA } from "$baseDir/subworkflows/local/metadata"

include { ARCHR_PROCESSING } from "$baseDir/subworkflows/local/archr_processing"
include { ARCHR_STAGE_PROCESSING } from "$baseDir/subworkflows/local/archr_stage_processing"

include { METADATA as METADATA_RNA } from "$baseDir/subworkflows/local/metadata"
include { ARCHR_INTEGRATE } from "$baseDir/subworkflows/local/archr_integration"

//
// SET CHANNELS
//

// set channel to reference folder containing fasta and gtf
Channel
    .value(params.reference)
    .set{ch_reference}


//
// WORKFLOW: Run main nf-core/downstream analysis pipeline
//
workflow NFCORE_DOWNSTREAM {

    METADATA( params.sample_sheet )

    // add gtf to cellranger output so can add annotations
    METADATA.out // METADATA.out: [[meta], [cellranger_output]]
        .combine(ch_reference)
        .map{[it[0], it[1] + it[2]]}
        .set {ch_metadata} // ch_metadata: [[meta], [cellranger_output, gtf]]

    // ARCHR: run processing + clustering + filtering + gene scores on full data
    ARCHR_PROCESSING ( ch_metadata ) //output = archr_filtered_full

    // ARCHR: run clustering + gene scores on individual stages
    ARCHR_STAGE_PROCESSING ( ARCHR_PROCESSING.out.output )

    // ATAC: combine stage data and full data
    ARCHR_STAGE_PROCESSING.out.output
        .combine(ARCHR_PROCESSING.out.output)
        .view()
        .set {ch_atac}

    // RNA: read in data
    //.set ch_rna
    METADATA_RNA( params.rna_sample_sheet )
        .set {ch_rna}
        //.view()
   
    // combine ATAC and RNA data
    ch_atac
        .mix(ch_rna)
        .groupTuple(by:0)
        .view()
        .set {ch_integrate}

    // ARCHR: Integrate
    ARCHR_INTEGRATE(ch_integrate)
}

//[[sample_id:ss8], [[rds_files], [ss8_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:ss4], [[rds_files], [ss4_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:HH5], [[rds_files], [HH5_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:HH7], [[rds_files], [HH7_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:HH6], [[rds_files], [HH6_cell_state_classification.RDS]], [[sample_id:FullData]], [FullData_Save-ArchR]]
//[[sample_id:full], [[seurat_label_transfer.RDS]]]

// /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/cc/b3b68fd1425e054dd6ae759ec13ea7/rds_files/HH6_Save-ArchR
// [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/cc/b3b68fd1425e054dd6ae759ec13ea7/rds_files/HH6_Save-ArchR, [sample_id:FullData], /camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/c4/133f89c7532f55fef898d7d089f2cb/rds_files/FullData_Save-ArchR]
// [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/cb/8d0fac51814fedee9d5febf5c38e49/rds_files/ss8_Save-ArchR, [[sample_id:FullData]], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/c4/133f89c7532f55fef898d7d089f2cb/rds_files/FullData_Save-ArchR]]
// [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/f7/bc54f69accf01c7f01966502879bf9/rds_files/ss4_Save-ArchR, [[sample_id:FullData]], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/c4/133f89c7532f55fef898d7d089f2cb/rds_files/FullData_Save-ArchR]]
// [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/67/c15bd31017b326b4d8c97277de573e/rds_files/HH5_Save-ArchR, [[sample_id:FullData]], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/c4/133f89c7532f55fef898d7d089f2cb/rds_files/FullData_Save-ArchR]]
// [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/ea/32136722dcfd8aa72f42d46b3bb19a/rds_files/HH7_Save-ArchR, [[sample_id:FullData]], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/c4/133f89c7532f55fef898d7d089f2cb/rds_files/FullData_Save-ArchR]]
// [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/cc/b3b68fd1425e054dd6ae759ec13ea7/rds_files/HH6_Save-ArchR, [[sample_id:FullData]], [/camp/svc/scratch/luscomben/hamrude/atac_neural_plate_border/NF-downstream_analysis/work/c4/133f89c7532f55fef898d7d089f2cb/rds_files/FullData_Save-ArchR]]
// [[sample_id:HH5], [[HH5_cell_state_classification.RDS]]]
// [[sample_id:HH6], [[HH6_cell_state_classification.RDS]]]
// [[sample_id:HH7], [[HH7_cell_state_classification.RDS]]]
// [[sample_id:ss4], [[ss4_cell_state_classification.RDS]]]
// [[sample_id:ss8], [[ss8_cell_state_classification.RDS]]]
// [[sample_id:FullData], [[seurat_label_transfer.RDS]]]

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
    NFCORE_DOWNSTREAM ()
}


/*
========================================================================================
    THE END
========================================================================================
*/
