process CUTADAPT {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    path(reads)
    path(mates)
    val(reads_adapter)
    val(mates_adapter)
    val(min_read_len)
    val(sample_name)

    output:
    path('*.converted_1.*'), emit: fq1
    path('*.converted_2.*'), optional: true, emit: fq2

    script:
    def adapter_outfile_info = "-a $reads_adapter -o TRIMMED.${reads.getBaseName()}"
    if (mates != "" && mates_adapter != ""){
        adapter_outfile_info += " -A $mates_adapter -p TRIMMED.${mates.getBaseName()}"
    }
    def input_reads = reads
    if (mates != ""){
        input_reads += " $mates"
    }

    """
    cutadapt -j 8 \\
    $adapter_outfile_info \\
    $input_reads
    """


    stub:
    """
    touch reads.fastq.gz
    """
}