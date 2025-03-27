process RMATS {
    label 'C2'
    container "xinglab/rmats:v4.2.0"

    input:
    path(gtf_annotation)
    path(sample_1_bams)
    val(read_length)
    val(read_type)
    val(strandedness)

    output:
    path('*.A3SS.*JC.txt'), emit: alternative_3_prime_splice_sites_jc
    path('*.A5SS.*JC.txt'), emit: alternative_5_prime_splice_sites_jc
    path('*.MXE.*JC.txt'), emit: mutually_exclusive_exons_jc
    path('*.RI.*JC.txt'), emit: retained_introns_jc
    path('*.SE.*JC.txt'), emit: skipped_exons_jc
    path('temp/*_read_outcomes_by_bam.txt'), emit: temp_read_outcomes
    path("$task.prefix/summary.txt"), emit: summary_file

    script:
    def ext_args = task.ext.args ?: ''
    """
    echo $sample_1_bams > sample_1.txt
    python /rmats/rmats.py \\
    --gtf $gtf_annotation \\
    --b1 sample_1.txt \\
    --od $task.prefix \\
    --tmp temp \\
    --nthread $task.cpus \\
    -t $read_type \\
    --libType $strandedness \\
    --readLength $read_length \\
    $ext_args \\
    && for i in ./$task.prefix/*.txt;
    do cp \$i ${task.prefix}.`basename \$i`; done
    """

}