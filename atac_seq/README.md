# ATAC-seq Analysis Pipeline

A production-grade, research-oriented ATAC-seq workflow for bulk ATAC-seq data with executable upstream/downstream scripts and reproducible environment specifications.

## Scientific Scope

This pipeline supports:

- Raw FASTQ quality control and adapter trimming.
- Alignment to a reference genome with duplicate handling and mitochondrial read filtering.
- Peak calling and consensus peak generation.
- Quantification across consensus peaks and DESeq2-based differential accessibility testing.
- Functional interpretation (motif enrichment, peak annotation, pathway context).

## Workflow Overview

```text
FASTQ
  └─> QC (FastQC/MultiQC)
      └─> Adapter trimming (Trim Galore/cutadapt)
          └─> Alignment (BWA-MEM2/Bowtie2)
              └─> Filtering + dedup + shift correction
                  └─> Peak calling (MACS3)
                      └─> Consensus peak set
                          └─> Count matrix
                              └─> Differential accessibility (DESeq2/edgeR)
                                  └─> Motif & functional analysis
```

## Directory Layout

- `upstream/`: raw-read processing through peak calling.
- `downstream/`: statistics, differential accessibility, and biological interpretation.

## Reproducibility & Quality Standards

- Fully pinned software environments (`environment.yml` or containerized execution).
- Mandatory metadata sheet with biological covariates and batch information.
- Multi-level QC checkpoints with explicit fail/warn thresholds.
- Separation of exploratory and confirmatory analyses.
- Versioned reference genome + blacklist files documented in run logs.

## How to Start

1. Prepare metadata from `metadata.template.tsv` (fill sample IDs and FASTQ paths).
2. Create environment from `environment.yml`.
3. Run `upstream/scripts/run_atac_upstream.sh` to generate BAM/peaks/count matrix.
4. Run `downstream/scripts/run_atac_downstream.sh` for DESeq2 analysis and QC outputs.

## Recommended Experimental Design

- Minimum biological replicates: **n >= 3 per condition**.
- For low-effect-size designs, target **n >= 4-6** per condition.
- Randomize library preparation and sequencing lanes to reduce confounding.

