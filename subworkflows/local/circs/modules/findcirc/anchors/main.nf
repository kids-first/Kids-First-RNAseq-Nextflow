process FINDCIRC_ANCHORS {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.auto_anchors.qfa'), emit: anchors

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /opt/circs_snake/scripts/pipelines/unmapped2anchors.py \\
    $reads \\
    > ${prefix}.auto_anchors.qfa
    """
}
