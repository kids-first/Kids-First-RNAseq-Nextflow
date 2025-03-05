process STAR_ALIGN {
    label 'process_star'
    container "pgc-images.sbgenomics.com/d3b-bixu/star:2.7.10a"

    input:
    path(genomeDir)
    val(readFilesCommand)
    val(read_groups)
    path(fastq_reads)
    val(pe_flag)
    val(outFileNamePrefix)

    output:
    path('*Log.progress.out'), emit: log_progress_out
    path('*Log.out'), emit: log_out
    path('*Log.final.out'), emit: log_final_out
    path('*Aligned.out.bam'), emit: genomic_bam_out
    path('*SJ.out.tab.gz'), emit: junctions_out
    path('*Aligned.toTranscriptome.out.bam'), emit: transcriptome_bam_out
    path('*Chimeric.out.sam'), optional: true, emit: chimeric_sam_out
    path('*Chimeric.out.junction'), emit: chimeric_junctions
    path('*ReadsPerGene.out.tab.gz'), emit: gene_counts

    script:
    def star_ext_args = task.ext.args ?: ''
    // Use pe_flag to help determine whether to pull two files at a time, or one for each RG
    def manifest_str = ""
    def i = 0
    
    read_groups.each{ v ->
            manifest_str +=  fastq_reads[i] + "\t" + (pe_flag ? fastq_reads[i+1] : "-" ) + "\t" + v + "\n";
            if (pe_flag){
                i += 2
            }
            else{
                i += 1
            }
        }
    """
    echo -e "$manifest_str" > star_reads_manifest.txt \\
    && tar -I pigz -xvf $genomeDir \\
    && STAR \\
    --genomeDir ./${genomeDir.getBaseName().replace(".tar", "")} \\
    --readFilesCommand $readFilesCommand \\
    --readFilesManifest star_reads_manifest.txt \\
    --outFileNamePrefix "${outFileNamePrefix}." \\
    --runThreadN $task.cpus \\
    $star_ext_args \\
    && pigz *tab
    """

}