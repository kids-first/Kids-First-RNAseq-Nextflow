process CIRCEXPLORER2_ANNOTATE {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circexplorer2:2.3.8"

    input:
    tuple val(meta), path(fusion_junctions)
    path(ref_fasta)
    path(refseq_annot)

    output:
    tuple val(meta), path('*known.txt'), emit: circs

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    CIRCexplorer2 annotate \\
    -b $fusion_junctions \\
    -g $ref_fasta \\
    -r $refseq_annot \\
    -o ${prefix}.circularRNA_known.txt \\
    $args
    """
}
