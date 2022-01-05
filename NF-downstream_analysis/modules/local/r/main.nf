

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process CELLRANGER_COUNT {
    tag "$meta.sample_id"
    label 'process_medium'

    //if (params.enable_conda) {
    //    exit 1, "Conda environments cannot be used when using the Cell Ranger tool. Please use docker or singularity containers."
    //}
    container "nfcore/cellranger:6.0.2"

    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), file('*')       , emit: outs
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
        Rscript ${params.script} --cores ${task.cpus} --runtype nextflow ${options.args}
        rm -r input
        """
}