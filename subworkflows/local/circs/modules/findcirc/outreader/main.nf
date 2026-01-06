process FINDCIRC_OUTREADER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(annotated_circs)

    output:
    path('*.findcirc.tsv'), emit: processed_circs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    perl /opt/circs_snake/scripts/f_c_outreader.pl \\
    --i $annotated_circs \\
    --o ${prefix}.findcirc.tsv \\
    $args 
    """
}
