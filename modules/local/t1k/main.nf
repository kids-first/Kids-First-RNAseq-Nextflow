process T1K {
    label 'M8'
    container 'pgc-images.sbgenomics.com/d3b-bixu/t1k:v1.0.5'

    input:
    tuple path(bam), path(bai)
    path(reference)
    path(gene_coordinates)

    output:
    path('*_aligned*.fa'), emit: aligned_fasta
    path('*_allele.tsv'), emit: allele_tsv
    path('*_allele.vcf'), emit: allele_vcf, optional: true
    path('*_candidate*.fq'), emit: candidate_fastqs
    path('*_genotype.tsv'), emit: genotype_tsv
    path('*_assign.tsv'), emit: read_assignments, optional: true

    script:
    def ext_args = task.ext.args ?: ''
    """
    run-t1k \\
    -b $bam \\
    -f $reference \\
    -c $gene_coordinates \\
    -o $task.ext.prefix \\
    $ext_args \\
    && sed -i '1s/^/gene_name\\tnum_diff_alleles\\tallele_1\\tabundance_1\\tquality_1\\tallele_2\\tabundance_2\\tquality_2\\tsecondary_alleles\\n/' *_genotype.tsv
    """

}