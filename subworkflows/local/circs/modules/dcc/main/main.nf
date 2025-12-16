process DCC_MAIN {
    label 'process_single'
    container "pgc-images.sbgenomics.com/danmiller/circs-dcc:0.1.0"

    input:
    tuple val(meta), path(paired_junctions), path(read1_junctions), path(read2_junctions)
    path(ref_fasta)
    path(refseq_bed)

    output:
    path('CircRNACount_clean'), emit: counts 
    path('CircCoordinates_clean'), emit: coordinates

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    /opt/scripts/pipelines/main.py \\
    $paired_junctions \\
    -mt1 $read1_junctions \\
    -mt2 $read2_junctions \\
    -an $refseq_bed \\
    -A $ref_fasta \\
    -T $task.cpus \\
    -D \\
    -fg \\
    -Pi \\
    -M \\
    -Nr 2 1 \\
    -N \\
    $args \\
    && sed '1d' CircRNACount > CircRNACount_clean \\\
    && sed '1d' CircCoordinates > CircCoordinates_clean
    """
}
