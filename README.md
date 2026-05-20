# CRISPR-Nanopore-QC-
A lightweight Streamlit web app for rapid quality control and sgRNA representation scoring of CRISPR screen populations using Oxford Nanopore sequencing.

# CRISPR-Nanopore-QC 🧬

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**CRISPR-Nanopore-QC** is a lightweight, fast web application designed to provide a rapid sanity check for CRISPR screen populations using Oxford Nanopore Technologies (ONT) sequencing. 

Waiting for deep Illumina sequencing to confirm library representation can cost valuable time and money. This tool allows researchers to sequence a PCR amplicon of their sgRNA library on a MinION/Flongle (yielding a small 5-10 MB FASTQ file) and immediately assess guide distribution, representation, and potential jackpotting before committing to next-generation sequencing.

## ✨ Features

* **Speed & Efficiency:** Optimized for free-tier hosting limits using rapid C-based Levenshtein distance calculations (`edlib`) instead of heavy brute-force alignments.
* **Indel-Tolerant Extraction:** Extracts the ~20bp sgRNA barcode by fuzzy-matching highly conserved 5' (e.g., U6) and 3' (e.g., scaffold) flanking sequences, elegantly handling Nanopore's characteristic indel profile.
* **Bidirectional Handling:** Automatically checks both the forward and reverse-complement strands to account for random pore transit orientations.
* **Quality Control:** Filters out junk reads using a customizable Phred quality score threshold.
* **Composite Ranking:** Ranks genes based on a composite score that balances total read depth with the fraction of expected guides successfully detected.

## 🚀 Live Demo
(https://b3yrslgizptbgufb9bersq.streamlit.app/#rapid-crispr-screen-qc-nanopore)

## 📋 Data Requirements

To run the analysis, you need two files:

### 1. Nanopore FASTQ File
A small FASTQ file (5-10 MB is typically sufficient, representing ~10,000 to 30,000 reads) from a PCR amplification of the sgRNA region. 

### 2. sgRNA Library Key (CSV)
A comma-separated values (`.csv`) file containing your library key. **The column headers must match exactly as shown below, with no trailing spaces.** 

| Guide_Seq | Gene | Guide_Num |
| :--- | :--- | :--- |
| AGCTACCATGCAGTACGTAG | SCGB1A1 | 1 |
| TCGATGCATCGATCGTACGA | SCGB1A1 | 2 |
| CCATGCAGTACGTAGCTAGC | CD4 | 1 |
| CGTAGCTAGCTAGCTAGCTA | CD4 | 2 |

* **`Guide_Seq`**: The DNA sequence of the guide (uppercase A, T, C, G). Do not include the PAM sequence.
* **`Gene`**: The target gene symbol or control identifier (e.g., `NonTargeting_1`).
* **`Guide_Num`**: An integer representing the guide's index for that specific gene.

## 💻 Local Installation

If you prefer to run the app locally rather than on the cloud:

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/yourusername/CRISPR-Nanopore-QC.git](https://github.com/yourusername/CRISPR-Nanopore-QC.git)
   cd CRISPR-Nanopore-QC

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
3. **Run the Streamlit app:**
   ```bash
   streamlit run app.py

How the Scoring Works
Due to PCR biases and sequencing noise, raw read counts alone can be misleading (e.g., a single "jackpot" guide dominating the reads for a gene). To provide a more robust representation metric, genes are ranked using a composite score:

Score = Log2(1 + Total Reads) * (Guides Detected / Expected Guides)

This formula explicitly rewards genes where multiple independent guides are detected, rather than genes with a single highly amplified sequence.

🛠️ Built With
Streamlit - The web framework

Biopython - FASTQ parsing and strand manipulation

edlib - Lightweight, super-fast C/C++ string matching

Plotly - Interactive data visualization
