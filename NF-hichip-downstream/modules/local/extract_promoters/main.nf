process EXTRACT_PROMOTERS {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple path(gtf_file), path(index_file)

    output:
    path "$output"       , emit: bed

    script:
    output = gtf_file.toString() - ".gtf" + ".bed"
    """
    extract_promoters.sh $gtf_file $index_file $output
    rm -r $gtf_file $index_file 
    """
}