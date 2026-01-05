process CIRCEXPLORER_MAIN {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(fusion_junctions)
    path(ref_fasta)
    path(refseq_annot)

    output:
    tuple val(meta), path('*circ.txt'), emit: circs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    python /opt/circs_snake/scripts/pipelines/CIRCexplorer.py \\
    -j $fusion_junctions \\
    -g $ref_fasta \\
    -r $refseq_annot \\
    -o $prefix \\
    $args
    """
}
