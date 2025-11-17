#!/bin/bash
set -euo pipefail

# 定义变量
MANIFEST_PATH="../manifest"
METADATA_PATH="../metadata.tsv"
CLASSIFIER_PATH="/public/home/mengpf/Project/16s/database/classifier_silva_13_8_99_V3-V4.qza"
THREADS=30
DATE=$(date +"%Y%m%d")
WORKDIR=$(pwd)

# 定义日志函数
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# 检查命令是否执行成功
check_success() {
    if [ $? -ne 0 ]; then
        log "Error: $1 failed"
        exit 1
    fi
}

# 检查必要文件
check_files() {
    for file in "$MANIFEST_PATH" "$METADATA_PATH" "$CLASSIFIER_PATH"; do
        if [ ! -f "$file" ]; then
            log "Error: Required file not found: $file"
            exit 1
        fi
    done
}

# 创建目录
create_directory() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
        log "Created directory: $1"
    fi
}

# 1. 导入数据
import_data() {
    log "Starting data import"
    create_directory "1.import_data"
    cd 1.import_data

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "$MANIFEST_PATH" \
        --output-path paired_end_demux_${DATE}.qza \
        --input-format PairedEndFastqManifestPhred33V2
    check_success "Data import"

    qiime demux summarize \
        --i-data paired_end_demux_${DATE}.qza \
        --o-visualization demux_summary_${DATE}.qzv
    check_success "Demux summarize"

    cd ..
    log "Data import completed"
}

# 2. DADA2 处理
run_dada2() {
    log "Starting DADA2 processing"
    create_directory "2.dada2"
    cd 2.dada2

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ../1.import_data/paired_end_demux_${DATE}.qza \
        --p-trim-left-f 19 \
        --p-trim-left-r 20 \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --o-table feature_table_${DATE}.qza \
        --o-representative-sequences rep_seqs_${DATE}.qza \
        --o-denoising-stats denoising_stats_${DATE}.qza \
        --p-n-threads $THREADS

           qiime feature-table summarize \
      --i-table feature_table_${DATE}.qza    \
      --o-visualization table.qzv
      
    qiime metadata tabulate \
        --m-input-file denoising_stats_${DATE}.qza \
        --o-visualization denoising-stats.qzv
        
    check_success "DADA2 denoise-paired"

    cd ..
    log "DADA2 processing completed"
}

# 3. 物种注释
run_taxonomy() {
    log "Starting taxonomy classification"
    create_directory "3.taxonomy"

    qiime feature-classifier classify-sklearn \
        --i-classifier "$CLASSIFIER_PATH" \
        --i-reads 2.dada2/rep_seqs_${DATE}.qza \
        --o-classification 3.taxonomy/taxonomy_${DATE}.qza
    check_success "Taxonomy classification"

    log "Taxonomy classification completed"
}

# 4. 导出数据
export_data() {
    log "Exporting data"
    create_directory "4.exports"

    # 导出ASV表
    qiime tools export \
        --input-path 2.dada2/feature_table_${DATE}.qza \
        --output-path 4.exports/feature_table
    check_success "Feature table export"

    # 导出物种注释
    qiime tools export \
        --input-path 3.taxonomy/taxonomy_${DATE}.qza \
        --output-path 4.exports/taxonomy
    check_success "Taxonomy export"

    # 导出代表序列
    qiime tools export \
        --input-path 2.dada2/rep_seqs_${DATE}.qza \
        --output-path 4.exports/rep_seqs
    check_success "Representative sequences export"

    # 转换为TSV格式
    biom convert -i 4.exports/feature_table/feature-table.biom \
                -o 4.exports/feature_table/feature-table.tsv --to-tsv
    check_success "BIOM to TSV conversion"

    log "Data export completed"
}

# 主函数
main() {
    log "Starting QIIME2 analysis pipeline"
    
    check_files
    import_data
    run_dada2
    run_taxonomy
    export_data
    
    log "Analysis pipeline completed successfully"
    log "Results can be found in 4.exports directory"
}

# 运行主函数
main