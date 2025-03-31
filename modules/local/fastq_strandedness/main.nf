process FASTQ_STRANDEDNESS {
    label 'C4'
    container "pgc-images.sbgenomics.com/d3b-bixu/stranded:1.1.0"

    input:
    tuple val(meta_reads), path(reads)
    path(annotation_gtf)
    path(kallisto_idx)
    val(nreads)

    output:
    path('fastq.strandness'), emit: result
    env(TOP_LEN), emit: top_read_len

    script:
    def reads_in = reads.size() == 2
        ? "--reads_1 ${reads[0]} --reads_2 ${reads[1]}"
        : "--reads_1 ${reads}"
    def pseudo_bam = "stranded_test_${reads[0].getBaseName() - ~/.f\w*q$/}/kallisto_strand_test/pseudoalignments.bam"
    """
    check_strandedness \\
    --gtf $annotation_gtf \\
    --kallisto_index $kallisto_idx \\
    $reads_in \\
    --nreads $nreads > fastq.strandness
    
    TOP_LEN=`samtools view $pseudo_bam \\
    | cut -f 10 \\
    | awk '{print length}' | sort | uniq -c | sort -nr -k1,1 | head -n 1`

    """

    stub:
    """
    touch fastq.strandness
    """
}