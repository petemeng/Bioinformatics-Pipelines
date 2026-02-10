# ATAC-seq Analysis Pipeline

A production-grade, bulk ATAC-seq workflow that emphasizes **reproducibility, traceable QC, and publication-ready statistics**.

## What Is Implemented

This pipeline is already executable end-to-end (not just a template):

- Upstream: FASTQ QC, alignment, filtering, deduplication, blacklist removal, MACS3 peak calling, consensus peak generation, featureCounts matrix.
- Downstream: DESeq2 differential accessibility, PCA QC plot, statistical summary, and session information export.

See:
- `upstream/scripts/run_atac_upstream.sh`
- `downstream/scripts/run_atac_downstream.sh`

## Quick Start (5 steps)

1. Create environment

```bash
conda env create -f environment.yml
conda activate atac-seq-pipeline
```

2. Prepare metadata

```bash
cp metadata.template.tsv metadata.tsv
# edit metadata.tsv: fill sample_id/group/replicate/batch/fastq_r1/fastq_r2
```

3. Run upstream

```bash
bash upstream/scripts/run_atac_upstream.sh \
  --metadata metadata.tsv \
  --reference /path/to/genome.fa \
  --blacklist /path/to/blacklist.bed \
  --outdir results_upstream \
  --threads 16 \
  --mapq 30 \
  --genome-size hs
```

4. Run downstream

```bash
bash downstream/scripts/run_atac_downstream.sh \
  --counts results_upstream/matrix/consensus_peak_counts.txt \
  --metadata metadata.tsv \
  --outdir results_downstream \
  --contrast TRT,CTRL
```

5. Check key outputs

- `results_upstream/consensus/consensus_peaks.bed`
- `results_upstream/matrix/consensus_peak_counts.txt`
- `results_downstream/differential/deseq2_results.tsv`
- `results_downstream/qc/pca_vst.png`

## Workflow Overview

```text
FASTQ
  -> FastQC / MultiQC
  -> bwa-mem2 alignment
  -> samtools filter/dedup (+ remove chrM)
  -> bedtools blacklist filtering
  -> MACS3 peaks per sample
  -> consensus peak merge
  -> featureCounts matrix (SAF)
  -> DESeq2 differential accessibility + PCA QC
```

## Project Structure

- `upstream/` : FASTQ -> BAM -> Peaks -> Matrix
- `downstream/` : Matrix -> Differential accessibility + QC
- `metadata.template.tsv` : metadata format reference
- `environment.yml` : pinned dependencies

## Scientific and QC Recommendations

- Biological replicates: ideally `>= 3` per group.
- Always inspect duplication rate, mtDNA fraction, and peak-level QC before interpretation.
- If strong batch effects are expected, extend model design to include batch covariates.
- Report both effect size (`log2FC`) and FDR, not only p-values.

## Current Scope and Limitations

- Scope: **bulk ATAC-seq** (paired-end) with group-level differential analysis.
- Current default model: `design = ~ group` (downstream script).
- Motif enrichment / peak-to-gene annotation are recommended as next extensions.

## Reproducibility Notes

- Scripts write run summaries and R session information.
- Reference genome and blacklist should be versioned and recorded in project logs.
- Use `--dry-run` first to validate file paths and arguments.
