process DCC_MAIN {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/danmiller/circtools:2.0.4"

    input:
    tuple val(meta), path(paired_junctions), path(read1_junctions), path(read2_junctions), path(sj_tabs)
    path(ref_fasta)
    path(refseq_bed)

    output:
    tuple val(meta), path('CircRNACount_clean'), emit: counts
    tuple val(meta), path('CircCoordinates_clean'), emit: coordinates

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    gzip -d *.SJ.out.tab.gz \\
    && echo $paired_junctions > samplesheet.txt \\
    && echo $read1_junctions > mate1.txt \\
    && echo $read2_junctions > mate2.txt \\
    && circtools detect \\
    @samplesheet.txt \\
    -mt1 @mate1.txt \\
    -mt2 @mate2.txt \\
    -an $refseq_bed \\
    -A $ref_fasta \\
    -T $task.cpus \\
    $args \\
    && sed '1d' CircRNACount > CircRNACount_clean \\
    && sed '1d' CircCoordinates > CircCoordinates_clean
    """
}
