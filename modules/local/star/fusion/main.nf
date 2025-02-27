process STAR_FUSION {
    label 'process_star'
    container "pgc-images.sbgenomics.com/d3b-bixu/star:fusion-1.10.1"

    input:
    path(genome_tar)
    path(Chimeric_junction)
    val(output_basename)

    output:
    path('*fusion_predictions.abridged.coding_effect.tsv'), emit: abridged_coding
    path('*.gz'), emit: chimeric_junction_compressed

    script:
    def star_ext_args = task.ext.args ?: ''
    """
    tar -I pigz -xvf $genome_tar \\
    && /usr/local/STAR-Fusion/STAR-Fusion \\
    --output_dir STAR-Fusion_outdir \\
    -J $Chimeric_junction \\
    $star_ext_args \\
    && mv STAR-Fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv ${output_basename}.STAR-1.10.1.fusion_predictions.abridged.coding_effect.tsv \\
    && pigz -c $Chimeric_junction > ${Chimeric_junction.getName()}.gz
    """

}