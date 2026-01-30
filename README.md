# Probe Enrichment Platform Analysis

This repository contains the processing pipeline for DNA probe libraries enriched by hybridizing against a specific target. The workflow begins with raw FASTQ files and utilizes a mapping strategy, where the probe motifs are count. 


---

## Getting Started

### Prerequisites
The following tools are required for the analysis:
* **SeqKit (v2.10.0):** For motif mapping and target identification.
* **FastQC (v0.11.9):** For quality control assessment.
* **Python/R:** (Insert versions used) for downstream statistical analysis.

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

### 3. PolyT Barcode Mapping
Unique barcode motifs are utilized for mapping to circumvent polyT non-homogeneity:
`seqkit locate --pattern-file barcodes.fa sample_R1.fastq.gz`
