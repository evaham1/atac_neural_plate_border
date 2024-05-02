process PYTHON {
    //tag "$meta.sample_id"
    label 'process_medium'

    container "alexthiery/seacells:latest"

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
    ${params.script} --input input/rds_files/ --ncores ${task.cpus} ${args} 
    rm -r input
    """
}