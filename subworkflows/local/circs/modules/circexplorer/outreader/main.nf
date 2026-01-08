process CIRCEXPLORER_OUTREADER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(circs)

    output:
    tuple val(meta), path('*cx.tsv'), emit: processed_circs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = meta.sample_name ?: "${meta.id}"
    """
    perl /opt/circs_snake/scripts/circexplorer1_out_reader.pl \\
    $circs \\
    ${prefix}.cx.tsv \\
    $sample_name
    """
}
