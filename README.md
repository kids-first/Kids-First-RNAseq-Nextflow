# Kids First RNAseq Nextflow
This repo is currently a dev project of converting our CWL production workflow to Nextflow.
Once development has complete, it will become a prod product.

<p align="center">
  <img src="docs/kids_first_logo.svg" alt="Kids First repository logo" width="660px" />
</p>
<p align="center">
  <a href="https://github.com/kids-first/Kids-First-RNAseq-Nextflow/blob/main/LICENSE"><img src="https://img.shields.io/github/license/kids-first/Kids-First-RNAseq-Nextflow.svg?style=for-the-badge"></a>
</p>

## Current State
The workflow is in alpha production shape.
Currently, you __cannot__ mix single and and paried end data with this workflow.
In the `test_inputs` dir, there are examples of various input situations that have been tested.

## Preprocess Reads Subworkflow:
The workflow takes in a mix of alignment files (BAM/CRAM) and fastq (single or paired end) and does the following:
### Alignment input
1. Split by read group
1. Create STAR read group strings
1. Convert to FASTQ
1. Run cutadapt if a cutadapt-related param is given
1. Return an object with STAR RG strings and related fastqs for downstream processing

### FASTQ input
1. Reformat to match object created by alignment input
1. Run cutadapt if a cutadapt-related param is given
### Result
Return an object with STAR RG strings and related fastqs coming from the alignment input and/or fastq input for downstream processing and the added_metadata with the following:
   - Paired end flag
   - Read length median
   - Read length std dev
   - Strandedness

## Align Analyze RNAseq
 - STAR Align
 - STAR Fusion
 - Arriba Fusion
 - T1K
 - RNASeQC
 - Kallisto

## annoFuse Subworkflow
 - Format Arriba
 - Annotate Arriba
 - Collate, filter, and annotate Arriba + Fusion results (annoFuse)

## T1K Subworkflow
 - Run T1K
 - Filter results
# DAGS

## Main WF
```mermaid
flowchart TB
    subgraph " "
    subgraph params
    v6["read_length_median"]
    v16["sample_id"]
    v48["hla_rna_gene_coords"]
    v8["read_length_stddev"]
    v20["reference"]
    v2["input_fastq_reads"]
    v32["readFilesCommand"]
    v34["FusionGenome"]
    v38["assembly"]
    v24["output_basename"]
    v22["reference_index"]
    v18["threads"]
    v4["is_paired_end"]
    v44["RNAseQC_GTF"]
    v12["max_reads"]
    v14["line_filter"]
    v10["strandedness"]
    v30["genomeDir"]
    v46["hla_rna_ref_seqs"]
    v26["gtf_anno"]
    v40["samtools_threads"]
    v36["fusion_annotator_tar"]
    v28["kallisto_idx"]
    v42["RSEM_genome"]
    v0["input_alignment_reads"]
    end
    v50([preprocess_reads])
    v56([align_analyze_rnaseq])
    v57([annofuse_subworkflow])
    v58([rmats_subworkflow])
    v0 --> v50
    v2 --> v50
    v4 --> v50
    v6 --> v50
    v8 --> v50
    v10 --> v50
    v12 --> v50
    v14 --> v50
    v16 --> v50
    v18 --> v50
    v20 --> v50
    v26 --> v50
    v28 --> v50
    v32 --> v56
    v34 --> v56
    v38 --> v56
    v40 --> v56
    v42 --> v56
    v44 --> v56
    v46 --> v56
    v16 --> v56
    v48 --> v56
    v50 --> v56
    v20 --> v56
    v22 --> v56
    v24 --> v56
    v26 --> v56
    v28 --> v56
    v30 --> v56
    v16 --> v57
    v36 --> v57
    v56 --> v57
    v24 --> v57
    v50 --> v58
    v56 --> v58
    v24 --> v58
    v26 --> v58
    end
```

## Preprocess Reads

```mermaid
flowchart TB
    subgraph preprocess_reads
    subgraph take
    v4["read_length_median"]
    v8["sample_id"]
    v9["threads"]
    v3["is_paired_end"]
    v11["annotation_gtf"]
    v5["read_length_stddev"]
    v10["reference"]
    v7["max_reads"]
    v1["input_fastq_reads"]
    v2["line_filter"]
    v6["strandedness"]
    v12["kallisto_idx"]
    v0["input_alignment_reads"]
    end
    v14([SAMTOOLS_HEAD_RG_CT])
    v18([SAMTOOLS_SPLIT])
    v20([SAMTOOLS_HEAD])
    v22([ALIGNMENT_PAIREDNESS])
    v25([SAMTOOLS_FASTQ])
    v30([FASTQ_STRANDEDNESS])
    v38([CUTADAPT])
    subgraph emit
    v42["added_metadata"]
    v41["cutadapt_stats"]
    v40["fastq_to_align"]
    end
    v0 --> v14
    v9 --> v18
    v10 --> v18
    v14 --> v18
    v2 --> v20
    v18 --> v20
    v14 --> v20
    v0 --> v22
    v7 --> v22
    v9 --> v22
    v10 --> v22
    v20 --> v25
    v22 --> v25
    v8 --> v25
    v9 --> v25
    v10 --> v25
    v1 --> v30
    v7 --> v30
    v25 --> v30
    v11 --> v30
    v12 --> v30
    v1 --> v38
    v25 --> v38
    v38 --> v40
    v38 --> v41
    v1 --> v42
    v5 --> v42
    v6 --> v42
    v25 --> v42
    v30 --> v42
    end
```

## Align Analyze RNAseq

```mermaid
flowchart TB
    subgraph align_analyze_rnaseq
    subgraph take
    v4["genome_tar"]
    v7["reference_index"]
    v12["sample_id"]
    v2["outFileNamePrefix"]
    v10["RSEMgenome"]
    v15["hla_rna_gene_coords"]
    v13["RNAseQC_GTF"]
    v3["input_fastq_reads"]
    v0["genomeDir"]
    v14["hla_rna_ref_seqs"]
    v1["readFilesCommand"]
    v8["gtf_anno"]
    v5["samtools_threads"]
    v16["added_metadata"]
    v9["assembly"]
    v6["reference_fasta"]
    v11["kallisto_idx"]
    end
    v23([STAR_ALIGN])
    v24([STAR_FUSION])
    v25([SAMTOOLS_SORT])
    v27([KALLISTO])
    v29([ARRIBA_FUSION])
    v30([ARRIBA_DRAW])
    v31([RSEM])
    v32([RNASEQC])
    v33([TAR_GZ])
    v34([T1K])
    v36([SAMTOOLS_VIEW])
    subgraph emit
    v47["RSEM_gene"]
    v49["RNASeQC_counts"]
    v46["RSEM_isoform"]
    v40["STAR_gene_count"]
    v50["kallisto_Abundance"]
    v41["STAR_junctions_out"]
    v51["t1k_genotype_tsv"]
    v37["genomic_bam_out"]
    v48["RNASeQC_Metrics"]
    v44["arriba_fusion_results"]
    v42["STAR_final_log"]
    v43["STARFusion_results"]
    v38["STAR_sorted_genomic_cram"]
    v45["arriba_fusion_viz"]
    v39["STAR_chimeric_junctions"]
    end
    v0 --> v23
    v16 --> v23
    v1 --> v23
    v2 --> v23
    v3 --> v23
    v2 --> v24
    v4 --> v24
    v23 --> v24
    v2 --> v25
    v5 --> v25
    v23 --> v25
    v16 --> v27
    v3 --> v27
    v11 --> v27
    v12 --> v27
    v16 --> v29
    v2 --> v29
    v6 --> v29
    v8 --> v29
    v9 --> v29
    v25 --> v29
    v2 --> v30
    v8 --> v30
    v9 --> v30
    v25 --> v30
    v29 --> v30
    v16 --> v31
    v2 --> v31
    v23 --> v31
    v10 --> v31
    v16 --> v32
    v25 --> v32
    v13 --> v32
    v32 --> v33
    v2 --> v33
    v2 --> v34
    v25 --> v34
    v14 --> v34
    v15 --> v34
    v6 --> v36
    v7 --> v36
    v25 --> v36
    v23 --> v37
    v36 --> v38
    v24 --> v39
    v23 --> v40
    v23 --> v41
    v23 --> v42
    v24 --> v43
    v29 --> v44
    v30 --> v45
    v31 --> v46
    v31 --> v47
    v32 --> v48
    v33 --> v49
    v27 --> v50
    v34 --> v51
    end
```

## annoFuse Subworkflow

```mermaid
flowchart TB
    subgraph annofuse_subworkflow
    subgraph take
    v3["rsem_expr_file"]
    v5["output_basename"]
    v2["fusion_annotator_tar"]
    v4["star_fusion_output_file"]
    v0["arriba_output_file"]
    v1["sample_name"]
    end
    v6([FORMAT_ARRIBA])
    v7([ANNOTATE_ARRIBA])
    v8([ANNOFUSE])
    subgraph emit
    v9["annofuse_filtered_fusions_tsv"]
    end
    v0 --> v6
    v5 --> v6
    v2 --> v7
    v5 --> v7
    v6 --> v7
    v1 --> v8
    v3 --> v8
    v4 --> v8
    v5 --> v8
    v7 --> v8
    v8 --> v9
    end
```

## rMATS Subworkflow

```mermaid
flowchart TB
    subgraph rmats_subworkflow
    subgraph take
    v1["sample_1_bams"]
    v3["read_type"]
    v4["strandedness"]
    v5["output_basename"]
    v0["gtf_annotation"]
    v2["read_length"]
    end
    v6([RMATS])
    v8([AWK_JC_FILTER])
    subgraph emit
    v9["rmats_filtered_jc"]
    end
    v0 --> v6
    v1 --> v6
    v2 --> v6
    v3 --> v6
    v4 --> v6
    v5 --> v6
    v6 --> v8
    v8 --> v9
    end
```