process R {
    //tag "$meta.sample_id"
    label 'process_medium'

    //if (params.enable_conda) {
    //    exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    //}
    container "rocker/rstudio:latest"

    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), file('*')       , emit: outs
    //path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export HDF5_USE_FILE_LOCKING=FALSE
    Rscript ${params.script} --cores ${task.cpus} --runtype nextflow ${args}
    rm -r input
    """
}
