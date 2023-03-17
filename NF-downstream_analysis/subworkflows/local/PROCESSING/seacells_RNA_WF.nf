#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Python scripts to run SEACells computation
include {R as SEURAT_TO_ANNDATA} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )

include {PYTHON as CALCULATE_SEACELLS} from "$moduleDir/../../../modules/local/python/main"               addParams(script: file("$moduleDir/../../../bin/seacells/SEACells_computation.py", checkIfExists: true) )
include {PYTHON as EXPORT_DATA_FROM_SEACELLS} from "$moduleDir/../../../modules/local/python/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/Export_data_from_AnnData.py", checkIfExists: true) )

include {R as META_TO_SEURAT} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seacells_meta_to_seurat.R", checkIfExists: true) )
include {R as SEURAT_TO_ANNDATA_PROCESSED} from "$moduleDir/../../../modules/local/r/main"               addParams(script: file("$moduleDir/../../../bin/data_conversion/seurat_to_h5ad.R", checkIfExists: true) )

//include {R as CHECK_SEACELLS} from "$baseDir/modules/local/r/main"               addParams(script: file("$baseDir/bin/seacells/Check_seacell_purity.R", checkIfExists: true) )


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate metacells

workflow SEACELLS_RNA_WF {
    take:
    input //[[sample_id:TransferLabels], [Processing/TransferLabels/3_peak_call/rds_files/TransferLabels_Save-ArchR]]

    main:

    input.set { ch_seurat }

    // Convert seurat to Anndata object
    SEURAT_TO_ANNDATA( input )

    //////// Run SEACells /////////
    CALCULATE_SEACELLS( SEURAT_TO_ANNDATA.out ) // Python script to calculate seacells on AnnData object
    EXPORT_DATA_FROM_SEACELLS( CALCULATE_SEACELLS.out ) //Python script to export data from Anndata objects as .csv

    // Process resulting metacells - need to input original seurat object and the anndata exported data
    ch_combined = EXPORT_DATA_FROM_SEACELLS.out
            .concat( ch_seurat )
            .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
            .collect()
            .map { [[sample_id:'Input'], it] } // [[meta], [rds1, rds2, rds3, ...]]
    META_TO_SEURAT( ch_combined ) // this outputs 2 seurat objects, one full object with metacells added and one summarised seurat
    
    // need to filter to only the condensed seurat object to convert into anndata 

    //SEURAT_TO_ANNDATA_PROCESSED( ch )

    emit:
    anndata = CALCULATE_SEACELLS.out
    exports = EXPORT_DATA_FROM_SEACELLS.out
    seurat = META_TO_SEURAT.out
    //processed_anndata = SEURAT_TO_ANNDATA_PROCESSED.out

}
