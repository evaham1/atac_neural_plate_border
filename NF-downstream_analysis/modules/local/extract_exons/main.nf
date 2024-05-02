process EXTRACT_EXONS {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path(gtf_file)

    output:
    path "$output"       , emit: bed

    script:
    output = gtf_file.toString() - ".gtf" + ".bed"
    """
    extract_exons.sh $gtf_file $output
    rm -r $gtf_file
    """
}