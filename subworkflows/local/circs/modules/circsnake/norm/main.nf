process CIRCSNAKE_NORM {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-snake:0.1.0"

    input:
    tuple val(meta), path(voted_circs), path(reads_per_sample)

    output:
    tuple val(meta), path('*.normed_voted.csv'), emit: normed_voted_circs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript /opt/circsnake/scripts/norm_a_voted_circs_df.R \\
    $voted_circs \\
    $reads_per_sample \\
    ${prefix}.normed_voted.csv
    """
}
