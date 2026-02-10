# ATAC-seq Upstream Pipeline (FASTQ -> Peaks -> Count Matrix)

## Objectives

Generate production-ready ATAC-seq artifacts:

- Filtered and deduplicated BAM files
- Per-sample MACS3 peak calls
- Project-level consensus peak set
- Consensus peak count matrix (featureCounts)
- MultiQC summary and run metadata logs

## Input Files

1. `metadata.tsv` with columns:
   - `sample_id`, `group`, `replicate`, `batch`, `fastq_r1`, `fastq_r2`
2. Reference genome FASTA (`genome.fa`) compatible with `bwa-mem2`
3. ENCODE blacklist BED

## Implemented Processing Steps

- FastQC for each paired FASTQ
- Alignment with `bwa-mem2 mem`
- BAM sorting/indexing (`samtools`)
- MAPQ + flag filtering, remove mitochondrial reads (`chrM`)
- Duplicate removal (`samtools markdup -r`)
- Blacklist filtering (`bedtools intersect -v`)
- Peak calling (`macs3 callpeak --call-summits`)
- Peak merge to consensus (`bedtools merge`)
- Count matrix generation (`featureCounts`, SAF annotation)
- MultiQC aggregation

## Run Command

```bash
bash scripts/run_atac_upstream.sh \
  --metadata metadata.tsv \
  --reference /path/to/genome.fa \
  --blacklist /path/to/blacklist.bed \
  --outdir results_upstream \
  --threads 16 \
  --mapq 30 \
  --genome-size hs
```

## Dry Run (validate paths/plan)

```bash
bash scripts/run_atac_upstream.sh \
  --metadata metadata.tsv \
  --reference /path/to/genome.fa \
  --blacklist /path/to/blacklist.bed \
  --outdir results_upstream \
  --dry-run
```

## Output Structure

- `results_upstream/qc/` : FastQC + MultiQC
- `results_upstream/bam/` : intermediate and final BAMs
- `results_upstream/peaks/` : sample MACS3 peaks
- `results_upstream/consensus/consensus_peaks.bed`
- `results_upstream/matrix/consensus_peak_counts.txt`
- `results_upstream/logs/` : run summary

## Quality Notes

- 建议最少每组 3 个生物学重复。
- `MAPQ>=30`、`FDR<0.05` 是常见起点，不应替代特定课题优化。
- `chrM` 过滤策略需与课题（核ATAC/全细胞）一致。

## Common Pitfalls

- Metadata sample IDs must uniquely map to FASTQ pairs.
- Ensure reference FASTA used for alignment matches the blacklist genome build.
- For large cohorts, consider per-sample parallel execution in a scheduler (SLURM/SGE).
