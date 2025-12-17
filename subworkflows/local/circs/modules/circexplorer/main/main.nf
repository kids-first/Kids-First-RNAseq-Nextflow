process CIRCEXPLORER_MAIN {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-cx:0.1.0"

    input:
    tuple val(meta), path(fusion_junctions)
    path(ref_fasta)
    path(refseq_annot)

    output:
    path('CIRCexplorer_circ.txt'), emit: circs

    script:
    def args = task.ext.args ?: ''
    """
    /opt/scripts/pipelines/CIRCexplorer.py \\
    -j $fusion_junctions \\
    -g $ref_fasta \\
    -r $refseq_annot \\
    $args
    """
}
