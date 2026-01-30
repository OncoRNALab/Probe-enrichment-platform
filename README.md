Probe Enrichment Platform Analysis
This repository contains the processing pipeline for DNA libraries enriched for KRAS and BRAF targets. The workflow begins with raw FASTQ files and utilizes a synthetic barcode mapping strategy to overcome sequencing challenges associated with low-diversity polyT tracts.

Project Overview
The core of this platform is a custom barcoding strategy. To prevent base-calling errors and registration failures caused by monochromatic polyT sequences on the Illumina NextSeq 2000, we replace low-diversity signals with high-diversity synthetic barcodes (Annex 4).

Getting Started
Prerequisites

The following tools are required for the analysis:

SeqKit (v2.10.0): For motif mapping and target identification.

FastQC (v0.11.9): For quality control assessment.

Python/R: (Add your version here) for downstream statistical analysis.

Data Preparation

The pipeline expects paired-end FASTQ files generated from an Illumina NextSeq 2000 (P1 kit).

Loading Concentration: 650 pM

PhiX Spike-in: 15%

Demultiplexing: Performed via Bcl2fastq v2.20.
