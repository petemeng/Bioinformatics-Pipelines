# 16S rRNA Amplicon Sequencing Analysis Pipeline

This pipeline provides a complete workflow for 16S rRNA gene amplicon sequencing data, from raw reads to final analysis and visualization.

## 🧬 Pipeline Overview

The analysis is divided into two main stages:

1.  **[Upstream Analysis](./upstream/README.md)**: This stage focuses on processing the raw sequencing reads. It includes steps like quality control (QC), primer trimming, merging paired-end reads, chimera removal, and generating an Amplicon Sequence Variant (ASV) or Operational Taxonomic Unit (OTU) table.
    
2.  **[Downstream Analysis](./downstream/README.md)**: This stage focuses on the statistical analysis and biological interpretation of the ASV/OTU table. It includes alpha diversity, beta diversity, taxonomic composition analysis, and differential abundance testing.

##  workflow Diagram

```
[Raw FASTQ Files]
       │
       ▼
┌───────────────────┐
│  Upstream Analysis  │
│ (DADA2 / QIIME2)  │
└───────────────────┘
       │
       ▼
[ASV/OTU Table + Taxonomy]
       │
       ▼
┌────────────────────┐
│ Downstream Analysis│
│   (Phyloseq in R)  │
└────────────────────┘
       │
       ▼
[Plots and Statistical Results]
```

## 🚀 How to Run

To execute this pipeline, please follow the stages in order:

1.  **Start with the [Upstream Analysis](./upstream/README.md)** to process your raw data.
2.  Once you have the output (ASV/OTU table), proceed to the **[Downstream Analysis](./downstream/README.md)** for interpretation.

## 🧪 Test Data

Sample data for this pipeline is located in the `./data` directory. This allows you to run the pipeline and verify that it is working correctly.
