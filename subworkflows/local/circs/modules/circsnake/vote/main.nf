process CIRCSNAKE_VOTE {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-parent:0.1.0"

    input:
    tuple val(meta), path(fc_matrix), path(cx_matrix), path(dcc_matrix)

    output:
    tuple val(meta), path('*circex_approved_by_all_three.csv'), emit: circex_voted
    tuple val(meta), path('*dcc_approved_by_all_three.csv'), emit: dcc_voted
    tuple val(meta), path('*find_circ_approved_by_all_three.csv'), emit: find_circ_voted
    tuple val(meta), path('*.tiff'), emit: venn

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    auto_voting_usable_hpc_outputs.R \\
    $fc_matrix \\
    $cx_matrix \\
    $dcc_matrix \\
    $prefix.
    """
}
