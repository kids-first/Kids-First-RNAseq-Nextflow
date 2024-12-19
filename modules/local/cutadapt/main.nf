process CUTADAPT {
    label 'process_medium'
    container "pgc-images.sbgenomics.com/d3b-bixu/cutadapt:3.4"

    input:
    path(reads)
    path(mates)
    val(reads_adapter)
    val(mates_adapter)
    val(min_read_len)

    output:
    path("TRIMMED.${reads.getName()}"), emit: fq1
    path("TRIMMED.${mates.getName()}"), optional: true, emit: fq2
    path("cutadapt_results.txt"), emit: cutadapt_metrics

    script:
    def adapter_outfile_info = "-a $reads_adapter -o TRIMMED.${reads.getName()}"
    if (mates != ".nextflow.log" && mates_adapter != ""){
        adapter_outfile_info += " -A $mates_adapter -p TRIMMED.${mates.getName()}"
    }
    def input_reads = reads
    if (mates.getName() != ".nextflow.log"){
        input_reads += " $mates"
    }

    """
    cutadapt -j 8 \\
    $adapter_outfile_info \\
    $input_reads \\
    > cutadapt_results.txt
    """


    stub:
    """
    touch TRIMMED.reads.fastq.gz
    """
}

workflow cutadapt {
    take:
    reads
    mates
    reads_adapter
    mates_adapter
    min_read_len
    
    main:
    CUTADAPT(reads, mates, reads_adapter, mates_adapter, min_read_len)
    emit:
    fq1 = CUTADAPT.out.fq1
    fq2 = CUTADAPT.out.fq2
    cutadapt_metrics = CUTADAPT.out.cutadapt_metrics
}

workflow  {
    reads = Channel.fromPath(params.input_fastq_reads)
    mates = Channel.fromPath(params.input_fastq_mates)
    reads_adapter = Channel.value(params.reads_adapter)
    mates_adapter = Channel.value(params.mates_adapter)
    min_read_len = Channel.value(params.cutadapt_min_read_len)
    cutadapt(reads, mates, reads_adapter, mates_adapter, min_read_len)
}