process CUTADAPT {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    // Single end or paired end, 1 or 2 files, no other possibilities for this module
    tuple val(meta_reads), path(reads, arity: '1..2')

    output:
    tuple val(meta_reads), path("TRIMMED.*"), emit: fastq_out
    tuple val(meta_reads), path("cutadapt_stats.txt"), emit: cutadapt_metrics

    script:
    def args = task.ext.args ?: ''
    """
    cutadapt -j $task.cpus \\
    $args \\
    $reads \\
    > cutadapt_stats.txt
    """


    stub:
    """
    touch TRIMMED.reads.fastq.gz
    """
}
