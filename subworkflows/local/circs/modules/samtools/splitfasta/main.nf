process SAMTOOLS_SPLITFASTA {
    label 'C8'
    container "pgc-images.sbgenomics.com/danmiller/samtools_parallel:1.20"

    input:
    tuple path(ref_fasta), path(ref_index)

    output:
    path('split_fasta'), emit: split_fasta

    script:
    """
    mkdir split_fasta
    cut -f1 $ref_index | parallel -j $task.cpus 'samtools faidx $ref_fasta {} > split_fasta/{}.fa'
    """
}
