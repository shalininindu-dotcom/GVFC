# Genomic Variant Frequency Calculator

A small Streamlit web application for calculating Variant Allele Frequency (VAF).

## Formula

VAF (%) = Variant Alleles / Total Alleles × 100

## Features

- Calculate the frequency of a single genomic variant
- Upload a CSV containing multiple variants
- Validate allele sequences using Biopython
- Display results in a table
- Download calculated results as CSV

## CSV format

The uploaded CSV must contain:

- `variant_id`
- `variant_alleles`
- `total_alleles`

Example:

```csv
variant_id,variant_alleles,total_alleles
rs123456,10,100
rs789012,25,200
```

## Run locally

Install the requirements:

```bash
python -m pip install -r requirements.txt
```

Then start the app:

```bash
python -m streamlit run app.py
```

If `python` is not recognized on Windows, try:

```bash
py -m pip install -r requirements.txt
py -m streamlit run app.py
```
