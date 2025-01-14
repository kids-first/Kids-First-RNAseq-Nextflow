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
The workflow is partially "finished" with the initial work of preprocessing input reads being completed.
Currently, you __cannot__ mix single and and pairied end data with this workflow.
In the `test_inputs` dir, there are examples of various input situations that have been tested.
Currently this includes:
 - [Paired end BAM with multiple RGs and cutadapt params](test_inputs/main_bam_in.json)
 - [Paired end BAM with multiple RGs, no cutadapt params](test_inputs/main_bam_inno.json)
 - [Paired end FASTQ with multiple RGs, and cutadapt params](test_inputs/main_fastq_in.json)
 - [Paired end FASTQ with multiple RGs, no cutadapt params](test_inputs/main_fastq_inno.json)
 - [Single end FASTQ with multiple RGs, and cutadapt params](test_inputs/main_fastq_single.json)
 - [Two paired end BAMs, one with two RGs, the other with one RG and cutadapt params](test_inputs/main_mixed_rg_ct_bam.json)
 - [One single end BAM, with one RG and cutadapt params](test_inputs/main_se_bam.json)

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
Return an object with STAR RG strings and related fastqs coming from the alignment input and/or fastq input for downstream processing