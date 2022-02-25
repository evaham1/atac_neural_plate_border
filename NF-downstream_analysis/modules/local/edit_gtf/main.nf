process EDIT_GTF {
    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), file('*')       , emit: outs

    script:
    """
    cat ./input/galgal6/tag_chroms.gtf | sed -e '/^#/! s/^/chr/' > temp.gtf
    rm -r input
    """
}
