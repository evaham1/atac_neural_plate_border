#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// python script to integrate seacells
include {PYTHON as INTEGRATE_SEACELLS} from "$baseDir/modules/local/python/main"               addParams(script: file("$baseDir/bin/Integration/SEACells_integration.py", checkIfExists: true) )

// r script to process integration maps that are outputted and checks on seurat object
include {R as PROCESS_INTEGRATION_OUTPUTS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/Integration/process_SEACell_integration_outputs.R", checkIfExists: true) )


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

workflow SEACELLS_INTEGRATING {
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

        INTEGRATE_SEACELLS.out.view()
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

        ch_labels_combined.view()

        PROCESS_INTEGRATION_OUTPUTS( ch_labels_combined )

    //emit integrated outputs:
    emit:
    raw_integration_output = INTEGRATE_SEACELLS.out
    processed_integration_output = PROCESS_INTEGRATION_OUTPUTS.out
}

// [[sample_id:HH6], [[plots, rds_files], csv_files, integrated_output]]
// [[sample_id:HH7], [[plots, rds_files], csv_files, integrated_output]]
// [[sample_id:HH5], [[a3/73e2f29913e4d9ae0738d711104b23/plots, a3/73e2f29913e4d9ae0738d711104b23/rds_files], b5/52af06cc3559a5d45a976102bc509f/csv_files, 6b/5063c856a50fa6a02c4e7ecd2a8946/integrated_output]]
// [[sample_id:ss8], [[7a/300bb19288c200cb37814a2b7b390c/plots, 7a/300bb19288c200cb37814a2b7b390c/rds_files], a8/f7a307efa1093759a27fd61ea09350/csv_files, 82/ebddad91b6894dc83d95a7958eec9c/integrated_output]]
// [[sample_id:ss4], [[f6/da89740bd86e6fb6f71e323043d6be/plots, f6/da89740bd86e6fb6f71e323043d6be/rds_files], 71/b91188feb85928e3c2096811627513/csv_files, ee/518f3f3b2279062e158b863392c05f/integrated_output]]

