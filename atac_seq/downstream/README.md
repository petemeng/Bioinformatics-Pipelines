# ATAC-seq Downstream Pipeline (Count Matrix -> Differential Accessibility)

## Objectives

对共识峰计数矩阵进行统计建模，输出可发表级别的差异可及性结果与样本QC图。

## Required Inputs

- `consensus_peak_counts.txt`（来自上游 `featureCounts`）
- `metadata.tsv`（至少含 `sample_id`, `group`）

## Implemented Analysis (DESeq2)

- 读取并解析 featureCounts 输出
- 按 `sample_id` 与 metadata 对齐样本
- 低表达峰过滤（至少2个样本 count>=10）
- `design = ~ group` 差异分析
- 输出 `deseq2_results.tsv`
- 输出 VST-PCA 图 (`qc/pca_vst.png`)
- 输出 `summary.txt` 与 `sessionInfo.txt`

## Run Command

```bash
bash scripts/run_atac_downstream.sh \
  --counts results_upstream/matrix/consensus_peak_counts.txt \
  --metadata metadata.tsv \
  --outdir results_downstream
```

指定比较方向（实验组,对照组）：

```bash
bash scripts/run_atac_downstream.sh \
  --counts results_upstream/matrix/consensus_peak_counts.txt \
  --metadata metadata.tsv \
  --outdir results_downstream \
  --contrast TRT,CTRL
```

## Dry Run

```bash
bash scripts/run_atac_downstream.sh \
  --counts results_upstream/matrix/consensus_peak_counts.txt \
  --metadata metadata.tsv \
  --outdir results_downstream \
  --dry-run
```

## Output Files

- `results_downstream/differential/deseq2_results.tsv`
- `results_downstream/qc/pca_vst.png`
- `results_downstream/logs/summary.txt`
- `results_downstream/logs/sessionInfo.txt`

## Statistical Notes

- 若存在批次效应，建议扩展设计公式到 `~ batch + group`。
- 建议报告 FDR、log2FC 及生物学效应解释，而不仅是 p 值。
