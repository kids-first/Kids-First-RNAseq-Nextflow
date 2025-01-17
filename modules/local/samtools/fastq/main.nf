process SAMTOOLS_FASTQ {
    label 'process_medium'
    container "staphb/samtools:1.20"

    input:
    tuple val(meta), path(input_align)
    path(cram_reference)
    val(threads)
    val(is_paired_end)

    output:
    tuple val(meta), path('*.converted*')

    script:
    def cram_ref_param =
        cram_reference != "" ? "--reference $cram_reference"
        : ''
    def pe_output_str = "-1 ${input_align.getBaseName()}.converted_1.fastq.gz -2 ${input_align.getBaseName()}.converted_2.fastq.gz -"
    def se_output_str = "-o ${input_align.getBaseName()}.converted.fastq.gz -"
    // is_paired_end is passed as an array to satisfy both single file and scatter
    def output_fastq =
        is_paired_end[0] ? pe_output_str
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
