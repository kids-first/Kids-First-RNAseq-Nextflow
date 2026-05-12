process FINDCIRC_OUTREADER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circexplorer2:2.3.8"

    input:
    tuple val(meta), path(annotated_circs)

    output:
    tuple val(meta), path('*.findcirc.tsv'), emit: processed_circs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    f_c_outreader.py \\
    --infile $annotated_circs \\
    --outfile ${prefix}.findcirc.tsv \\
    $args 
    """
}
