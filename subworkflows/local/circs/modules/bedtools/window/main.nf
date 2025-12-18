process BEDTOOLS_WINDOW {
    label 'process_single'
    container "staphb/bedtools:2.31.1"

    input:
    tuple val(meta), path(intervals1)
    path(intervals2)

    output:
    tuple val(meta), path('*.bed'), emit: windows

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    bedtools window \\
    -a $intervals1 \\
    -b $intervals2 \\
    -w 1 \\
    $args \\
    > ${prefix}.bed
    """
}
