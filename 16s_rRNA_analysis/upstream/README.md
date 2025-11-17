# 16S Pipeline: Upstream Analysis

This document provides detailed instructions for the upstream processing of 16S rRNA amplicon sequencing data. The goal is to convert raw FASTQ files into a feature table (ASV or OTU).

## 🛠️ Dependencies

Before running the script, ensure you have the following software installed and available in your PATH. Please add the specific tools your script uses.

*   **QIIME 2**: `(version X.X)`
*   **DADA2**: `(version X.X)`
*   **FastQC**: `(version X.X)`
*   **Cutadapt**: `(version X.X)`
*   *(Add any other tools required by your script)*

**Example Installation (using Conda):**
```bash
# Example for QIIME 2
# conda install ...
```

## 📜 Script

The primary script for this stage is:
*   `scripts/16s_upstream.sh`

## ▶️ Usage

1.  **Prepare your data**:
    Place your raw, paired-end `.fastq.gz` files in the `../data/` directory.

2.  **Modify the script (if necessary)**:
    Open `scripts/16s_upstream.sh` and adjust any variables at the top of the file, such as input directory, primer sequences, or truncation lengths.
    ```bash
    # --- User-configurable variables ---
    INPUT_DIR="path/to/your/data"
    OUTPUT_DIR="path/to/output"
    FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
    REVERSE_PRIMER="GACTACHVGGGTATCTAATCC"
    # ------------------------------------
    ```

3.  **Execute the script**:
    Run the script from the `16s_rRNA_analysis` directory:
    ```bash
    bash ./upstream/scripts/16s_upstream.sh
    ```

## 📤 Expected Output

Upon successful completion, the script will generate the following files in your specified output directory:

*   `otu_table.qza`: The main feature table.
*   `taxonomy.qza`: Taxonomic assignments for each feature.
*   `rep-seqs.qza`: Representative sequences for each feature.
*   *(Add any other important output files)*

These files will be the primary input for the **[Downstream Analysis](../downstream/README.md)**.
