process EDIT_VALIDPAIRS {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path(output_file)          , emit: validpairs

    script:
    output_file = input_file.toString() - ".allValidPairs" + "_edited.allValidPairs"
    """
    edit_validpairs.sh $input_file $output_file
    """
}

    

