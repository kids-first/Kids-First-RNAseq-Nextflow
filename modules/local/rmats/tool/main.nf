process RMATS {
    label 'process_annofuse'
    container "xinglab/rmats:v4.2.0"

    input:
    path(gtf_annotation)
    path(sample_1_bams)
    val(read_length)
    val(read_type)
    val(strandedness)
    val(output_basename)

    script:
    def ext_args = task.ext.args ?: ''
    """
    echo $sample_1_bams > sample_1.txt
    python /rmats/rmats.py \\
    --gtf $gtf_annotation \\
    --b1 sample_1.txt \\
    --od $output_basename \\
    --tmp temp \\
    --nthread 4 \\
    -t $read_type \\
    --libType $strandedness \\
    --readLength $read_length \\
    $ext_args
    """

}