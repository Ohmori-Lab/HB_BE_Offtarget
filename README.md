## Overall Description

This repository contains a set of bioinformatics tools designed to process and analyze genetic data, specifically focusing on the analysis of GUIDE-seq off-target regions, protein sequence logos, and DNA variant generation with protein translation. The tools cater to various stages of genetic analysis, from annotating genomic regions to generating protein sequence logos and simulating DNA sequence variants. Below is a summary of the three scripts included in this repository:

+ **DNA Variant Generator and Protein Sequence Translator:**
This tool allows for the creation of DNA sequence variants by applying specific substitutions (A->G, T->C, or both) in a given range of positions. The modified DNA sequences are then translated into their corresponding protein sequences using the standard codon table. The results are saved in a CSV file, detailing both the DNA variants and their corresponding protein sequences.

+ **Protein Sequence Logo Generator:**
This script generates a sequence logo from a set of protein sequences. Using a Position Specific Scoring Matrix (PSSM), it visualizes the conservation of amino acids across a collection of aligned protein sequences. The output is a graphical representation (sequence logo) saved in PNG format, highlighting amino acid frequency at each position in the sequence.

+ **GUIDE-seq Off-Target Region Annotator:**
This tool is designed for annotating off-target regions identified in GUIDE-seq experiments. It takes raw off-target data, normalizes genomic coordinates, and annotates the regions (exons, introns, intergenic regions) using a reference GTF annotation file (e.g., GENCODE v46 for hg38). The script generates detailed Excel summaries and visual plots to facilitate further analysis.

## Key Features:

+ **Comprehensive Genetic Data Processing:** Annotates off-target regions, generates sequence logos, and simulates DNA sequence variants with protein translation.

+ **Output Formats:** Generates both Excel and CSV files for downstream analysis, along with visual outputs like PNG plots for protein sequence logos.

+ **Flexible and Customizable:** Accepts user-defined input files, substitution types, and allows for the use of custom annotation files (e.g., GTF).

+ **Error Handling:** The scripts provide detailed error messages and validation checks to ensure input files are correctly formatted.

## Requirements:

+ **GUIDE-seq Off-Target Region Annotator:** Requires `pandas`, `numpy`, `matplotlib`, and `xlsxwriter` (optional for Excel output).

+ **Protein Sequence Logo Generator:** Requires `pandas`, `numpy`, `matplotlib`, and `logomaker`.

+ **DNA Variant Generator and Protein Sequence Translator:** Uses only Python's standard libraries.

## Applications:

+ **GUIDE-seq Analysis:** Annotating and visualizing off-target regions from CRISPR-based GUIDE-seq experiments.

+ **Sequence Logo Generation:** Analyzing and visualizing conserved regions in aligned protein sequences, useful in protein evolution and structure studies.

+ **Variant Simulation and Translation:** Simulating genetic variants and translating them to proteins, useful in genetic variant analysis and mutation studies.

These scripts provide a powerful toolkit for bioinformatics researchers working with genetic sequence data, from variant analysis to visualization of sequence conservation.
