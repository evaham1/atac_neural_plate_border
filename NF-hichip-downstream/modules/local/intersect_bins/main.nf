process INTERSECT_BINS {

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"
    
    input:
    tuple path(bins), path(features)

    output:
    path "$output"       , emit: intersected_bed

    script:
    output = gtf_file.toString() - ".bed" + "_bins_intersected.bed"
    """
    bedtools intersect -a $bins -b $features -wa -wb > $output
    """
}