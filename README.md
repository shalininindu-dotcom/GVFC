# Genomic Variant Frequency Calculator

Streamlit application for count-based variant frequency calculations and aligned FASTA sequence comparison.

## Features
- Single variant frequency calculator.
- CSV dataset frequency calculator.
- FASTA analysis using Biopython `Bio.SeqIO`.
- Reference sequence selection and position-by-position comparison.
- Alternate allele counts, frequencies, results table, and CSV download.

## FASTA assumptions
- Input must be a pre-aligned multiple-sequence FASTA file.
- All sequences must have equal length.
- Each FASTA record is treated as one observed sequence/allele.
- The selected record is the reference.
- A, C, G and T are counted; N and gaps are excluded from the denominator.
- Positions are reported using 1-based numbering.

The FASTA frequency is a sequence-level alternate allele proportion. It is not automatically diploid clinical VAF. This educational tool does not perform read mapping, variant calling from raw sequencing reads, or clinical interpretation.

## Run locally
```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```
On Windows, use `py -m` if `python` is not recognised.
