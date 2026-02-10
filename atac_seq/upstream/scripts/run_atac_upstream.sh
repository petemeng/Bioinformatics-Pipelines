#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage:
  $0 --metadata <metadata.tsv> --reference <genome.fa> --blacklist <blacklist.bed> \
     --outdir <outdir> [--threads 16] [--mapq 30] [--genome-size hs] [--dry-run]

Required metadata columns:
  sample_id    fastq_r1    fastq_r2    group    replicate

Description:
  Bulk ATAC-seq upstream workflow (FASTQ -> filtered BAM -> peaks -> count matrix).
  Steps include FastQC/MultiQC, alignment (bwa-mem2), filtering/dedup, MACS3,
  consensus peak merge, and featureCounts quantification.
USAGE
}

THREADS=8
MAPQ=30
GENOME_SIZE="hs"
DRY_RUN=0
METADATA=""
REFERENCE=""
BLACKLIST=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --metadata) METADATA="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --blacklist) BLACKLIST="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --genome-size) GENOME_SIZE="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

log() { echo "[$(date +'%F %T')] $*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }

run_cmd() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    log "[DRY-RUN] $*"
  else
    log "$*"
    eval "$@"
  fi
}

check_cmd() {
  command -v "$1" >/dev/null 2>&1 || fail "Required command not found: $1"
}

[[ -n "$METADATA" && -n "$REFERENCE" && -n "$BLACKLIST" && -n "$OUTDIR" ]] || {
  usage; fail "missing required arguments"
}
[[ -f "$METADATA" ]] || fail "metadata file not found: $METADATA"
[[ -f "$REFERENCE" ]] || fail "reference file not found: $REFERENCE"
[[ -f "$BLACKLIST" ]] || fail "blacklist file not found: $BLACKLIST"

if [[ "$DRY_RUN" -eq 0 ]]; then
  for c in fastqc multiqc bwa-mem2 samtools bedtools macs3 featureCounts awk; do
    check_cmd "$c"
  done
fi

mkdir -p "$OUTDIR"/{logs,qc/raw_qc,bam,tmp,peaks,consensus,matrix}

if ! head -n1 "$METADATA" | tr '\t' '\n' | grep -Fxq "sample_id"; then
  fail "metadata missing required header: sample_id"
fi
if ! head -n1 "$METADATA" | tr '\t' '\n' | grep -Fxq "fastq_r1"; then
  fail "metadata missing required header: fastq_r1"
fi
if ! head -n1 "$METADATA" | tr '\t' '\n' | grep -Fxq "fastq_r2"; then
  fail "metadata missing required header: fastq_r2"
fi

log "ATAC upstream started"
log "metadata=$METADATA reference=$REFERENCE blacklist=$BLACKLIST threads=$THREADS mapq=$MAPQ dry_run=$DRY_RUN"

# Skip header and iterate samples.
tail -n +2 "$METADATA" | while IFS=$'\t' read -r sample_id group replicate batch fastq_r1 fastq_r2 rest; do
  [[ -n "${sample_id:-}" ]] || continue
  [[ -f "$fastq_r1" ]] || fail "R1 FASTQ not found for $sample_id: $fastq_r1"
  [[ -f "$fastq_r2" ]] || fail "R2 FASTQ not found for $sample_id: $fastq_r2"

  bam_prefix="$OUTDIR/bam/${sample_id}"
  nodup_bam="${bam_prefix}.filtered.nodup.bam"

  run_cmd "fastqc -t $THREADS -o $OUTDIR/qc/raw_qc $fastq_r1 $fastq_r2"
  run_cmd "bwa-mem2 mem -t $THREADS $REFERENCE $fastq_r1 $fastq_r2 | samtools sort -@ $THREADS -o ${bam_prefix}.sorted.bam"
  run_cmd "samtools index ${bam_prefix}.sorted.bam"
  run_cmd "samtools view -@ $THREADS -b -q $MAPQ -F 1804 ${bam_prefix}.sorted.bam | samtools view -@ $THREADS -b -h - | awk 'BEGIN{OFS="\t"} /^@/ || \$3!="chrM"' | samtools view -@ $THREADS -b -o ${bam_prefix}.filtered.bam"
  run_cmd "samtools sort -@ $THREADS -o ${bam_prefix}.filtered.sorted.bam ${bam_prefix}.filtered.bam"
  run_cmd "samtools markdup -@ $THREADS -r ${bam_prefix}.filtered.sorted.bam $nodup_bam"
  run_cmd "samtools index $nodup_bam"

  # Remove blacklist regions and call peaks.
  run_cmd "bedtools intersect -v -abam $nodup_bam -b $BLACKLIST > ${bam_prefix}.clean.bam"
  run_cmd "samtools index ${bam_prefix}.clean.bam"
  run_cmd "macs3 callpeak -t ${bam_prefix}.clean.bam -f BAMPE -g $GENOME_SIZE -n ${sample_id} --outdir $OUTDIR/peaks --keep-dup all --call-summits"
done

run_cmd "multiqc -o $OUTDIR/qc $OUTDIR/qc/raw_qc $OUTDIR/peaks"
run_cmd "cat $OUTDIR/peaks/*_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > $OUTDIR/consensus/consensus_peaks.bed"
run_cmd "awk \'BEGIN{OFS=\"\t\"; print \"GeneID\tChr\tStart\tEnd\tStrand\"} {print \"peak_\"NR,\$1,\$2+1,\$3,\".\"}\' $OUTDIR/consensus/consensus_peaks.bed > $OUTDIR/consensus/consensus_peaks.saf"
run_cmd "featureCounts -a $OUTDIR/consensus/consensus_peaks.saf -F SAF -o $OUTDIR/matrix/consensus_peak_counts.txt -p -T $THREADS $OUTDIR/bam/*.clean.bam"

cat > "$OUTDIR/logs/run_summary.txt" <<SUMMARY
metadata=$METADATA
reference=$REFERENCE
blacklist=$BLACKLIST
threads=$THREADS
mapq=$MAPQ
genome_size=$GENOME_SIZE
dry_run=$DRY_RUN
SUMMARY

log "ATAC upstream completed: $OUTDIR"
