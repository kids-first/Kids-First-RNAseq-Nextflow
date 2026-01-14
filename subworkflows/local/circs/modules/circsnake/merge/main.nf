process CIRCSNAKE_MERGE {
    label 'process_single'
    container "python:3.13"

    input:
    tuple val(meta), path(cx_normed), path(dcc_normed), path(fc_normed)

    output:
    tuple val(meta), path('*.merged.csv'), emit: merged_circs 

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_normed_voted_circs.py \\
    $cx_normed \\
    $dcc_normed \\
    $fc_normed \\
    ${prefix}.normed_voted.merged.csv
    """
}
