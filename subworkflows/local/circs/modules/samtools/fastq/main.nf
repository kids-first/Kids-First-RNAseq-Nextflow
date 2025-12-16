process SAMTOOLS_FASTQ {
    label 'C8'
    container "pgc-images.sbgenomics.com/d3b-bixu/samtools:1.15.1"

    input:
    tuple val(meta), path(input_align)
    path(cram_reference)

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz")      , optional:true, emit: fastq
    tuple val(meta), path("*_singleton.fastq.gz")  , optional:true, emit: singleton

    script:
    def cram_ref_param = cram_reference != "" ? "--reference $cram_reference" : ''
    def se_output_str = "-s ${input_align.getBaseName()}.singletons.fastq.gz -1 ${input_align.getBaseName()}.converted_1.fastq.gz" 
    def pe_output_str = se_output_str + " -2 ${input_align.getBaseName()}.converted_2.fastq.gz"
    def output_fastq = meta.is_paired_end ? pe_output_str : se_output_str
    """
    samtools collate \\
    $input_align \\
    -O \\
    --output-fmt SAM \\
    -f \\
    -r 50000 \\
    -@ $task.cpus \\
    $cram_ref_param \\
    | samtools fastq \\
    -@ $task.cpus \\
    -O \\
    $output_fastq
    """

    stub:
    """
    touch reads.fastq.gz
    """
}
