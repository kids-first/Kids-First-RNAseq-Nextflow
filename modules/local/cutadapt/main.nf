process CUTADAPT {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    tuple val(meta_reads), path(reads)
    tuple val(meta_mates), path(mates)

    output:
    tuple val(meta_reads), path("TRIMMED.*"), emit: fastq_out
    tuple val(meta_reads), path("cutadapt_stats.txt"), emit: cutadapt_metrics

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
