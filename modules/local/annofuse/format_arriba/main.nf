process ANNOFUSE {
    label 'process_annofuse'
    container "pgc-images.sbgenomics.com/d3b-bixu/annofuse:0.92.0"

    input:
    path(arriba_fusion_file)
    val(sample_name)

    output:
    path('*arriba_formatted.tsv'), emit: formatted_fusion_tsv

    script:
    """
    formatArribaFusionCalls.R \\
    --fusionfile $arriba_fusion_file
    --tumorID $sample_name
    --outputfile ${sample_name}.arriba_formatted.tsv

    """

}