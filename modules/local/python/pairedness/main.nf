process ALIGNMENT_PAIREDNESS {
    label 'process_low'
    container "quay.io/biocontainers/pysam:0.22.0--py310h41dec4a_0"
    input:
    path input_reads
    path input_reference
    val max_reads
    val output_filename
    val cpus


    output:
    val identifier

    script:
    """"
    alignmentfile_pairedness.py
    """"
    stub:
    """
    
    """
}