process ALIGNMENT_PAIREDNESS {
    label 'process_low'
    container "quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0"

    input:
    path input_reads
    path input_reference
    val max_reads
    val output_filename
    val threads


    output:
    path "${output_filename}", emit: identifier

    script:
    """
    alignmentfile_pairedness.py \\
    --input_reads $input_reads \\
    --input_reference $input_reference \\
    --max_reads $max_reads \\
    --threads $threads \\
    > $output_filename
    """
    stub:
    """
    touch pairedness.txt
    """
}
