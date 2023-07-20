process HOMER_MOTIF_ENRICHMENT {

    container "rocker/conda_homer:0.1"

    input:
    tuple val(meta), path(input_file)
    # needs to be a bed file of peaks and the genome index file

    output:
    tuple val(meta), path(output_file)          , emit: validpairs
    # needs to be the folder name for outputs

    script:
    output_file = input_file.toString() - ".allValidPairs" + "_edited.allValidPairs"
    """
    findMotifsGenome.pl ERpeaks.txt hg18 ER_MotifOutput/ -size given
    """
}

    

