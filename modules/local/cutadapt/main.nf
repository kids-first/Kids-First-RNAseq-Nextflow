process CUTADAPT {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    path(reads)
    path(mates)

    output:
    path("TRIMMED.*"), emit: fastq_out
    path("cutadapt_stats.txt"), emit: cutadapt_metrics

    script:
    def args = task.ext.args ?: ''
    def read_files = mates ? reads + " " + mates : reads
    """
    cutadapt -j 8 \\
    $args \\
    $read_files \\
    > cutadapt_stats.txt
    """


    stub:
    """
    touch TRIMMED.reads.fastq.gz
    """
}
