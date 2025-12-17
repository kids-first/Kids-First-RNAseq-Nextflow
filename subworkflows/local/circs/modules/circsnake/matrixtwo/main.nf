process CIRCSNAKE_MATRIXTWO {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-snake:0.1.0"

    input:
    tuple val(meta), path(circ_matrix)
    path(micrornas)
    path(coding_circs)
    path(hallmarks)
    path(ensembl)

    output:
    tuple val(meta), path('*.mat2'), emit: matrix

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /opt/circsnake/scripts/matrixtwo_V4.pl \\
    --i $circ_matrix \\
    --o ${prefix}.mat2 \\
    --m $micrornas \\
    --c $coding_circs \\
    --h $hallmarks \\
    --e $ensembl \\
    --n /opt/circsnake/scripts/read_mapping.pl \\
    --excl_cb 1
    """
}
