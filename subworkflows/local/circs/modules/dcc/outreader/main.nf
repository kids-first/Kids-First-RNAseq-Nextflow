process DCC_OUTREADER {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(annotated_counts), path(cleaned_coordinates)

    output:
    path('CircRNACount_clean'), emit: counts 
    path('CircCoordinates_clean'), emit: coordinates

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = $meta.sample_name ?: "${meta.id}"
    """
    perl /opt/circs_snake/scripts/dcc_outreader.pl \\
    $annotated_counts \\
    $cleaned_coordinates \\
    ${prefix}.dcc.tsv \\
    $sample_name
    """
}
