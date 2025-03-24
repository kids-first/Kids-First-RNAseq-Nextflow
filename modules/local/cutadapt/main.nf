process CUTADAPT {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    tuple val(meta_reads), path(reads)
    val(output_basename)

    output:
    tuple val(meta_reads), path("TRIMMED.*"), emit: fastq_out
    tuple val(meta_reads), path("cutadapt_stats.txt"), emit: cutadapt_metrics

    script:
    def args = task.ext.args ?: ''
    """
    cutadapt -j $task.cpus \\
    $args \\
    $reads \\
    > ${output_basename}.cutadapt_stats.txt
    """


    stub:
    """
    touch TRIMMED.reads.fastq.gz
    """
}
