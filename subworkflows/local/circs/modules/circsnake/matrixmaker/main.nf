process CIRCSNAKE_MATRIXMAKER {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(found_circs)
    path(known_circs)
    path(refseq)

    output:
    tuple val(meta), path('*.mat1'), emit: matrix

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    """
    cat $found_circs > all_tsvs.tx

    perl /opt/circs_snake/scripts/matrixmaker-V4.pl \\
    --threads $task.cpus \\
    --infile all_tsvs.tx \\
    --circ_bed_file $known_circs \\
    --gene_mapping_file $refseq \\
    --outfile ${prefix}.mat1 \\
    $args
    """
}
