process FINDCIRC_ANCHORS {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-findcirc:0.1.0"

    input:
    tuple val(meta), path(reads)

    output:
    path('*.auto_anchors.qfa'), emit: anchors

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /opt/circ_snake/scripts/pipelines/unmapped2anchors.py \\
    $reads \\
    > ${prefix}.auto_anchors.qfa
    """
}
