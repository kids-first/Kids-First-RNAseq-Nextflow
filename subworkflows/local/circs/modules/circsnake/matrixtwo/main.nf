process CIRCSNAKE_MATRIXTWO {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

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
    def args = task.ext.args ?: ""
    """
    perl /opt/circs_snake/scripts/matrixtwo_V4.pl \\
    --infile $circ_matrix \\
    --mirRNA_file $micrornas \\
    --circbank_coding_file $coding_circs \\
    --hallmark_mapping_file $hallmarks \\
    --ensembl_file $ensembl \\
    --mapping_script /opt/circs_snake/scripts/read_mapping.pl \\
    --outfile ${prefix}.mat2 \\
    $args
    """
}
