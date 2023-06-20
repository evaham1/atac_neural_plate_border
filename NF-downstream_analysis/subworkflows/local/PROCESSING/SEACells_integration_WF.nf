#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// python script to integrate seacells
include {PYTHON as INTEGRATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/Integration/SEACells_integration.py", checkIfExists: true) )

// r script to create consensus label transfer map and checks thes labels on SEACells ATAC seurat object
include {R as LABEL_TRANSFER} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/SEACell_integration_label_transfer.R", checkIfExists: true) )


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow SEACELLS_INTEGRATING_WF {
    take:
    RNA_SEACells //SEACELLS_RNA_WF.out.seacells_anndata_processed_classified
    ATAC_SEACells //SEACELLS_ATAC_WF.out.seacells_anndata_processed_classified
    ATAC_SEACells_seurat //SEACELLS_ATAC_WF.out.seacells_seurat_processed_classified
    ATAC_SEACells_map //SEACELLS_ATAC_WF.out.seacell_outputs_named

    main:

        // Prep input channel

        //RNA_SEACells.view()
        // [[sample_id:HH6], [c8/0355d43ba613bc2159e10758c26572/plots, c8/0355d43ba613bc2159e10758c26572/rds_files]]
        // [[sample_id:HH5], [61/920e02828bc1d30d01745280da08cb/plots, 61/920e02828bc1d30d01745280da08cb/rds_files]]
        // [[sample_id:ss8], [64/bc0c8d3f045086356bc00a7ebde6ff/plots, 64/bc0c8d3f045086356bc00a7ebde6ff/rds_files]]
        // [[sample_id:ss4], [fa/1d64469d41d5a69f74d4d706bf5681/plots, fa/1d64469d41d5a69f74d4d706bf5681/rds_files]]
        // [[sample_id:HH7], [d0/6b148d2543b6b195b8ee2088fad0b7/plots, d0/6b148d2543b6b195b8ee2088fad0b7/rds_files]]

        RNA_SEACells
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()[0]] }
            .set {ch_RNA}

        //ch_RNA.view()
        // [[sample_id:HH6], c8/0355d43ba613bc2159e10758c26572/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:HH5], 61/920e02828bc1d30d01745280da08cb/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:HH7], d0/6b148d2543b6b195b8ee2088fad0b7/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:ss4], fa/1d64469d41d5a69f74d4d706bf5681/rds_files/AnnData_RNA.h5ad]
        // [[sample_id:ss8], 64/bc0c8d3f045086356bc00a7ebde6ff/rds_files/AnnData_RNA.h5ad]
        
        ATAC_SEACells
            .map{ meta, data -> [meta, data.findAll{it =~ /rds_files/}[0].listFiles()[0]] }
            .set {ch_ATAC}
        // ch_ATAC.view()
        // [[sample_id:HH6], 6d/4467b6b5de1ef8f33c2edc45d871b6/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:HH7], 4f/2aef89dec2263316c6a5645b259204/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:HH5], 78/c6ea5ab5e6e2cc3e0b5ccd6763f948/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:ss4], a4/9d8ea368614ff7badd70b2135e7046/rds_files/AnnData_ATAC.h5ad]
        // [[sample_id:ss8], 4a/a274c2d504b915e6a9fb8369584796/rds_files/AnnData_ATAC.h5ad]

        ch_RNA
            .concat( ch_ATAC )
            .groupTuple( by:0 )
            //.view() //[[sample_id:ss4], [fa/1d64469d41d5a69f74d4d706bf5681/rds_files/AnnData_RNA.h5ad, 4a/72711280da408db03e302a5926456a/rds_files/AnnData_ATAC.h5ad]]
            .set {ch_seacells_to_integrate}

        // Run integration

        INTEGRATE_SEACELLS( ch_seacells_to_integrate )

        // Process integration outputs

        //ATAC_SEACells_seurat.view()
        // [[sample_id:HH6], [plots, rds_files]]
        // [[sample_id:HH7], [plots, rds_files]]
        // [[sample_id:ss8], [plots, rds_files]]
        // [[sample_id:ss4], [plots, rds_files]]
        // [[sample_id:HH5], [plots, rds_files]]

        //ATAC_SEACells_map.view()
        // [[sample_id:HH7], csv_files]
        // [[sample_id:HH6], csv_files]
        // [[sample_id:ss4], csv_files]
        // [[sample_id:HH5], csv_files]
        // [[sample_id:ss8], csv_files]

        //INTEGRATE_SEACELLS.out.view()
        // [[sample_id:HH6], [integrated_output]]
        // [[sample_id:HH5], [integrated_output]]
        // [[sample_id:ss4], [integrated_output]]
        // [[sample_id:HH7], [integrated_output]]
        // [[sample_id:ss8], [integrated_output]]

        ch_labels_combined = ATAC_SEACells_seurat
            .concat( ATAC_SEACells_map )
            .concat( INTEGRATE_SEACELLS.out )
            .groupTuple( by:0 ) // all 3 outputs need to be grouped by stage
            .map{ meta, data -> [meta, data.flatten()] } // remove extra [] around outputs

        //ch_labels_combined.view()
        // [[sample_id:HH7], [plots, rds_files, csv_files, integrated_output]]
        // [[sample_id:HH6], [plots, rds_files, csv_files, integrated_output]]
        // [[sample_id:ss4], [plots, rds_files, csv_files, integrated_output]]
        // [[sample_id:HH5], [plots, rds_files, csv_files, integrated_output]]
        // [[sample_id:ss8], [plots, rds_files, csv_files, integrated_output]]

        LABEL_TRANSFER( ch_labels_combined )

    //emit integrated outputs:
    emit:
    raw_integration_output = INTEGRATE_SEACELLS.out
    processed_integration_output = LABEL_TRANSFER.out
}