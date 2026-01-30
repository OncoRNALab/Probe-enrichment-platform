# Probe Enrichment Platform Analysis

This repository contains the processing pipeline for DNA probe libraries enriched via hybridization against specific genomic targets. By starting directly from raw FASTQ files, the workflow quantifies the experimental "pulldown" through a motif-based mapping strategy. Instead of traditional alignment to a reference genome, this platform identifies and counts specific probe motifs within the sequencing reads using SeqKit v2.10.0.

## Getting Started

### Prerequisites
The following tools are required for the analysis:
* **SeqKit (v2.10.0):** For motif mapping and target identification.
* **R (2025.05.1+513):**  for downstream statistical analysis.

### Data Preparation
The pipeline expects paired-end FASTQ files generated from an Illumina NextSeq 2000 (P1 kit). 

| Parameter | Value |
| :--- | :--- |
| **Loading Concentration** | 650 pM |
| **PhiX Spike-in** | 15% |
| **Demultiplexing** | Bcl2fastq v2.20 |

---

## Analysis Workflow

### 1. Quality Control
Assess read quality using FastQC:
`fastqc *.fastq.gz -o ./qc_reports/`

### 2. Target Identification
Probe motifs for *KRAS* and *BRAF* (referenced in Annex 2 and 3) are mapped directly using `SeqKit`:
`seqkit locate --pattern-file motifs.fa sample_R1.fastq.gz`

Unique barcode motifs are utilized for mapping to circumvent polyT (referenced in Annex 4) non-homogeneity:
`seqkit locate --pattern-file barcodes.fa sample_R1.fastq.gz`

### 3. Intermediate output
A locate file that contains to location of the motifs in the reads

| seqID | patternName | pattern | strand | start | end | matched |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| VH00559...54550 | Probe_167 | TTGGCTA... | + | 1 | 17 | TTGGCTA... |

### 4. Enrichment Quantification
After the mapping is complete, a custom Bash script processes the locate results to generate final counts. 

To maximize specificity and account for the known architecture of the library, we only count motifs that align at the beginning of the read:

Target Probes (KRAS/BRAF): Only matches starting at position 1 and ending at 17 are counted.

PolyT Barcodes: Only matches within the 1 to 7 range are counted.

### 5. Probe counts

| Probe_ID | Count |
| :--- | :--- |
| **Probe_167** | 1450 |




