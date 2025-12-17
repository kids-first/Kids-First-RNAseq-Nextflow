process CIRCSNAKE_MATRIXMAKER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-snake:0.1.0"

    input:
    tuple val(meta), path(found_circs)
    path(known_circs)
    path(refseq)

    output:
    path('*.mat1'), emit: matrix

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /opt/circsnake/scripts/matrixmaker-V4.pl \\
    --i $found_circs \\
    --o ${prefix}.mat1 \\
    --c $known_circs \\
    --g $refseq
    """
}
