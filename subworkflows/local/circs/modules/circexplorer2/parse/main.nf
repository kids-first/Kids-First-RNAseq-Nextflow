process CIRCEXPLORER2_PARSE {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circexplorer2:2.3.8"

    input:
    tuple val(meta), path(chimeric_junctions)

    output:
    tuple val(meta), path('*.bed'), emit: fusion_junctions
    tuple val(meta), path('*.cx_star_parse.log'), emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    CIRCexplorer2 parse \\
    -t $meta.aligner \\
    -b ${prefix}.fusion_junctions.bed \\
    $args \\
    $chimeric_junctions \\
    > ${prefix}.cx_star_parse.log
    """
}
