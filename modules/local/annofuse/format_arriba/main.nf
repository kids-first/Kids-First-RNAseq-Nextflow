process FORMAT_ARRIBA {
    label 'process_annofuse'
    container "pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0"

    input:
    path(arriba_fusion_file)
    val(output_basename)

    output:
    path('*arriba_formatted.tsv'), emit: formatted_fusion_tsv

    script:
    """
    formatArribaFusionCalls.R \\
    --fusionfile $arriba_fusion_file \\
    --outputfile ${output_basename}.arriba_formatted.tsv
    """

}