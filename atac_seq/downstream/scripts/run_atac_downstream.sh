#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  $0 --counts <consensus_peak_counts.txt> --metadata <metadata.tsv> \
     --outdir <outdir> [--contrast TRT,CTRL] [--dry-run]

Description:
  Run downstream ATAC-seq differential accessibility analysis with DESeq2.
  Input count file should be featureCounts-style output from upstream script.
USAGE
}

COUNTS=""
METADATA=""
OUTDIR=""
CONTRAST=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --counts) COUNTS="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --contrast) CONTRAST="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

log() { echo "[$(date +'%F %T')] $*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }

[[ -n "$COUNTS" && -n "$METADATA" && -n "$OUTDIR" ]] || { usage; fail "missing required arguments"; }
[[ -f "$COUNTS" ]] || fail "count file not found: $COUNTS"
[[ -f "$METADATA" ]] || fail "metadata file not found: $METADATA"

mkdir -p "$OUTDIR"/{qc,differential,logs}

if [[ "$DRY_RUN" -eq 1 ]]; then
  log "[DRY-RUN] Would run DESeq2 analysis with counts=$COUNTS metadata=$METADATA contrast=$CONTRAST"
  exit 0
fi

command -v Rscript >/dev/null 2>&1 || fail "Rscript not found"

R_SCRIPT="$OUTDIR/logs/run_deseq2_atac.R"
cat > "$R_SCRIPT" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Need: counts metadata outdir [contrast]")
counts_file <- args[[1]]
metadata_file <- args[[2]]
outdir <- args[[3]]
contrast_arg <- ifelse(length(args) >= 4, args[[4]], "")

meta <- read.delim(metadata_file, check.names = FALSE)
required <- c("sample_id", "group")
if (!all(required %in% colnames(meta))) stop("metadata missing required columns: sample_id/group")

fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
# featureCounts format: first 6 columns are annotation, remaining are sample BAM paths.
if (ncol(fc) < 7) stop("featureCounts file appears invalid")
count_mat <- fc[, 7:ncol(fc)]
rownames(count_mat) <- fc[[1]]
colnames(count_mat) <- basename(gsub("\\.clean\\.bam$", "", colnames(count_mat)))

meta <- meta[match(colnames(count_mat), meta$sample_id), ]
if (any(is.na(meta$sample_id))) stop("sample_id mismatch between counts and metadata")
meta$group <- factor(meta$group)

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(count_mat)), colData = meta, design = ~ group)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds)

if (contrast_arg != "") {
  parts <- strsplit(contrast_arg, ",")[[1]]
  if (length(parts) != 2) stop("contrast must be: test,control")
  res <- results(dds, contrast = c("group", parts[[1]], parts[[2]]))
} else {
  lv <- levels(meta$group)
  if (length(lv) < 2) stop("Need at least 2 groups")
  res <- results(dds, contrast = c("group", lv[[2]], lv[[1]]))
}

res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, file.path(outdir, "differential", "deseq2_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)))
pca_df <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], group = meta$group)
p <- ggplot(pca_df, aes(PC1, PC2, color = group, label = Sample)) + geom_point(size = 3) + theme_bw()
ggsave(file.path(outdir, "qc", "pca_vst.png"), p, width = 6, height = 5)

sig <- subset(res_df, !is.na(padj) & padj < 0.05)
writeLines(c(
  paste0("total_peaks_tested=", nrow(res_df)),
  paste0("significant_peaks_fdr0.05=", nrow(sig))
), file.path(outdir, "logs", "summary.txt"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "logs", "sessionInfo.txt"))
RSCRIPT

log "Running DESeq2 downstream analysis"
Rscript "$R_SCRIPT" "$COUNTS" "$METADATA" "$OUTDIR" "$CONTRAST"
log "Downstream ATAC analysis completed: $OUTDIR"
