process HOMER_MOTIF_ENRICHMENT {

    container "rocker/conda_homer:0.1"

    input:
    tuple path(input_file), path(fasta_file)
    //needs to be a bed file of peaks and the genome fasta file

    output:
    path(output_file)          , emit: motif_enrichments
    //dont know what this is, probably doesnt matter as not going to pipe output to anything else

    script:
    """
    split_bed_run_homer_motif_enrichment.sh $input_file $fasta_file
    """
}

    

