# Probe Enrichment Platform Analysis

This repository contains the processing pipeline for DNA probe libraries enriched via hybridization against specific genomic targets. By starting directly from raw FASTQ files, the workflow quantifies the experimental "pulldown" through a motif-based mapping strategy. Instead of traditional alignment to a reference genome, this platform identifies and counts specific probe motifs within the sequencing reads using SeqKit v2.10.0.

A script example is available in Example/

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


### Computational Implementation (HPC)

The analysis was performed on a **High-Performance Computing (HPC) environment** to leverage parallel processing of multiple FASTQ files. 

* **Environment:** Linux-based HPC Cluster.
* **Scripting:** Modular **Bash scripting** was used to automate the transition from raw FASTQ files to final probe counts.
* **Job Scheduling:** (Optional: Add if you used Slurm/PBS, e.g., "Jobs were submitted via Slurm to handle batch processing.")
* **Scalability:** The pipeline is designed to process multiple library pools simultaneously, significantly reducing the wall-clock time for large-scale enrichment experiments.

---

## Analysis Workflow

### 1. Quality Control
Assess read quality using FastQC:
`fastqc *.fastq.gz -o ./qc_reports/`

### 2. Target Identification
Probe motifs for *KRAS* and *BRAF* (referenced in Input/Annex 2 and 3) are mapped directly using `SeqKit`:
`seqkit locate --pattern-file motifs.fa sample_R1.fastq.gz`

Unique barcode motifs are utilized for mapping to circumvent polyT (referenced in Input/Annex 4) non-homogeneity:
`seqkit locate --pattern-file barcodes.fa sample_R1.fastq.gz`

### 3. Intermediate output
A locate file (referenced in Output) that contains the location of the probe motifs or barcodes in the reads

| seqID | patternName | pattern | strand | start | end | matched |
| :--- | :--- | :--- | :---: | :---: | :---: | :--- |
| VH00559...54550 | Probe_167 | TTGGCTA... | + | 1 | 17 | TTGGCTA... |

### 4. Enrichment Quantification
After the mapping is complete, a custom Bash script processes the locate results to generate final counts (Referenced in Ouput). 

To maximize specificity and account for the known architecture of the library, we only count motifs that align at the beginning of the read:

Target Probes (KRAS/BRAF): Only matches starting at position 1 and ending at 17 are counted.

PolyT Barcodes: Only matches within the 1 to 7 range are counted.

### 5. Probe counts

| Probe_ID | Count |
| :--- | :--- |
| **Probe_167** | 1450 |


### 6. Downstream analysis

Each sample in the dataset has an associated **probe count file**.  
These files contain raw sequencing counts. The metadata of each sample can be coupled to these counts to generate a complete, annotated dataset.

A typical metadata file contains the following fields:

- Probe_ID  
- Salt concentration  
- Temperature  
- MT percentage  
- WT percentage  
- Probe sequence  
- Position of the substitution  
- Substitution type  
- Mismatch position  
- Highlight: the probe type  

Example of the structure of a final combined dataset:

| Probe_ID | Salt | Temp | MT | WT | Sequence | Substitution | Mismatch_pos | Raw_Count | Norm_MT | Norm_WT |
|---------|------|------|----|----|------------------|-------------|-------------|-----------|---------|---------|
| Probe_01 | 600 | 60°C | 100 | 0 | atcttgcctacgccaca | C→A | 17 | 515 | 0.507 | 2.780 |
| Probe_01 | 600 | 60°C | 50 | 50 | atcttgcctacgccaca | C→A | 17 | 683 | 0.433 | 6.777 |
| Probe_01 | 600 | 60°C | 25 | 75 | atcttgcctacgccaca | C→A | 17 | 396 | 0.377 |10.385 |

Each row corresponds to a single measurement of a specific probe under a defined combination of laboratory conditions.

The coupling of raw count files with their corresponding metadata can be performed using custom R scripts.

As you can see, the normalized counts can also be calculated from the raw counts.  
This normalization procedure is described in detail in the section **"Free-energy Penalty Calculation"** in the paper.

Final datasets for the main experiments (**KRAS, BRAF, and PolyT**) are provided in the `output/` directory.








