process ALIGN_TO_FASTQ {
    label 'process_medium'
    container "staphb/samtools:1.20"

    input:
    path(input_align)
    path(cram_reference)
    val(sample_id)
    val(threads)
    val(is_paired_end)

    output:
    path('*.converted_1.*'), emit: fq1
    path('*.converted_2.*'), optional: true, emit: fq2

    script:
    def cram_ref_param =
        cram_reference != "" ? "--reference $cram_reference"
        : ''
    def pe_output_str = "-1 ${sample_id}.converted_1.fastq.gz -2 ${sample_id}.converted_2.fastq.gz -"
    def se_output_str = "-o ${sample_id}.converted_1.fastq.gz -"
    def output_fastq =
        is_paired_end ? pe_output_str
        : se_output_str
    """
    samtools sort \\
    -m 1G \\
    -n \\
    -O SAM \\
    -@ $threads \\
    $cram_ref_param \\
    $input_align \\
    | samtools fastq \\
    -c 2 \\
    -@ $threads \\
    $output_fastq
    """


    stub:
    """
    touch reads.fastq.gz
    """
}

workflow align_to_fastq {
    take:
    input_align
    cram_reference
    sample_id
    threads
    is_paired_end
    
    main:
    ALIGN_TO_FASTQ(input_align, cram_reference, sample_id, threads, is_paired_end)
    emit:
    fq1 = ALIGN_TO_FASTQ.out.fq1
    fq2 = ALIGN_TO_FASTQ.out.fq2
}

workflow  {
    input_align = Channel.fromPath(params.unsorted_bam)
    cram_reference = Channel.fromPath(params.reference)
    sample_id = Channel.from(params.sample_id)
    threads = Channel.value(params.threads)
    is_paired_end = params.is_paired_end
    align_to_fastq(input_align, cram_reference, sample_id, threads, is_paired_end)
    
}