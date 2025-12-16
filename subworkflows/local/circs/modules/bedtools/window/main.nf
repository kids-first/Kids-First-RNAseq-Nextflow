process BEDTOOLS_WINDOW {
    label 'process_single'
    container "biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1"

    input:
    tuple val(meta), path(intervals1), path(intervals2)

    output:
    path('*.bed'), emit: windows

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
