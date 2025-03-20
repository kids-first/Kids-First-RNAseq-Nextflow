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
Currently, you __cannot__ mix single and and paired end data with this workflow.
In the `test_inputs` dir, there are examples of various input situations that have been tested.
It can now run on CAVATICA, you can push this app by:
 - [Install sbpack](https://docs.cavatica.org/docs/bring-nextflow-apps-to-cavatica#about-sbpack)
 - Using command `sbpack_nf --profile {your_profile} --appid {username}/{project}/kfdrc-rnaseq-nextflow --workflow-path /path/to/this/repo/Kids-First-RNAseq-Nextflow/ --entrypoint main.nf --sb-schema sb_nextflow_schema.yaml`

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

## [Preprocess Reads](docs/dags/PREPROCESS_READS_SUBWF.md)

## [Align Analyze RNAseq](docs/dags/ALIGN_ANALYZE_RNASEQ_SUBWF.md)

## [annoFuse Subworkflow](docs/dags/ANNOFUSE_SUBWORKFLOW.md)

## [rMATS Subworkflow](docs/dags/RMATS_SUBWORKFLOW.md)
