process CUTADAPT {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    tuple val(meta), path(reads)

    output:
    path("TRIMMED.*"), emit: fastq_out
    path("cutadapt_results.txt"), emit: cutadapt_metrics

    script:
    def args = task.ext.args ?: ''
    """
    cutadapt -j 8 \\
    $args \\
    $reads
    > cutadapt_results.txt
    """


    stub:
    """
    touch TRIMMED.reads.fastq.gz
    """
}
