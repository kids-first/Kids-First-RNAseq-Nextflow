process CIRCEXPLORER_STARPARSE {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(chimeric_junctions)

    output:
    tuple val(meta), path('*.fusion_junction.txt'), emit: fusion_junctions
    tuple val(meta), path('*.cx_star_parse.log'), emit: log

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /opt/circs_snake/scripts/pipelines/star_parse.py \\
    $chimeric_junctions \\
    ${prefix}.fusion_junction.txt \\
    > ${prefix}.cx_star_parse.log
    """
}
