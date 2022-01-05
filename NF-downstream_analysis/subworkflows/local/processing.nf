//
// Run preprocessing R scripts
//

include {R_SIGNAC as PREPROCESSING} from "$baseDir/modules/local/r_signac/main"               addParams(options: modules['PREPROCESSING'],
                                                                                                script: file("$baseDir/bin/1_preprocessing.R", checkIfExists: true) )


workflow PROCESSING {
    input = [
        [ id:'test'],
        "$baseDir/../output/NF-luslab_sc_multiomic/hh7_1_cellranger_atac/outs/*", checkIfExists: true)
    ]

    PREPROCESSING( input )

}
