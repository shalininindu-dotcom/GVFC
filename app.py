
# ============================================================
# GENOMIC VARIANT FREQUENCY CALCULATOR
# ============================================================
# PERSON A: Web interface and user interaction
# PERSON B: Data processing, Biopython analysis and calculations
# ============================================================


# ============================================================
# SHARED SECTION: IMPORTING LIBRARIES
# ============================================================

import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
from io import StringIO


# ============================================================
# PERSON B: DATA PROCESSING AND CALCULATION FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# PERSON B - FUNCTION 1: SINGLE VARIANT CALCULATION
# ------------------------------------------------------------

def calculate_single_frequency(variant_id, allele, variant_count, total): #total is total no.of alleles checked

    if not variant_id.strip(): #if varient id is null
        raise ValueError("Please enter a variant ID.")

    if variant_count > total: #making sure numbers make sense
        raise ValueError("Variant alleles cannot exceed total alleles.")

    allele = allele.strip().upper()

    if allele and any(base not in "ACGTN" for base in allele):
        raise ValueError("Allele should contain only A, C, G, T, or N.")

    # Create a Biopython sequence object
    if allele:
        Seq(allele)

    # Calculate variant allele frequency
    frequency = variant_count / total * 100

    # Return processed results
    return {
        "Variant ID": variant_id.strip(),
        "Variant Allele": allele or "Not specified",
        "Variant Alleles": int(variant_count),
        "Total Alleles": int(total),
        "Frequency (%)": round(frequency, 2)
    }


# ------------------------------------------------------------
# PERSON B - FUNCTION 2: CSV FILE PROCESSING
# ------------------------------------------------------------

def process_csv(csv_file):

    # Read uploaded CSV file using Pandas
    df = pd.read_csv(csv_file) #df = dataframe

    # Define required CSV columns
    required = {
        "variant_id",
        "variant_alleles",
        "total_alleles"
    }

    # Check for missing columns
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Missing column(s): " + ", ".join(sorted(missing))
        ) #stops the code and lists the missing column names cleanly using

    # Convert allele counts into numeric values
    df["variant_alleles"] = pd.to_numeric(df["variant_alleles"], errors="coerce")

    df["total_alleles"] = pd.to_numeric(df["total_alleles"], errors="coerce")

    # Check for missing or non-numeric values
    if df[["variant_alleles", "total_alleles"]].isna().any().any(): #marks every blank cell as true and 
        raise ValueError("Allele counts must be numeric.")

    # Check total allele counts
    if (df["total_alleles"] <= 0).any():
        raise ValueError("Total allele counts must be greater than zero.")

    # Check negative variant counts
    if (df["variant_alleles"] < 0).any():
        raise ValueError("Variant counts cannot be negative.")

    # Check whether variant count exceeds total count
    if (df["variant_alleles"] > df["total_alleles"]).any():
        raise ValueError("Variant counts cannot exceed total counts.")

    # Calculate frequency for every row
    df["frequency_percent"] = (
        df["variant_alleles"] / df["total_alleles"] * 100
    ).round(2)

    return df


# ------------------------------------------------------------
# PERSON B - FUNCTION 3: FASTA FILE PARSING
# ------------------------------------------------------------

def parse_fasta(fasta_bytes):

    # Read FASTA records using Biopython SeqIO
    fasta_text = fasta_bytes.decode("utf-8-sig")
    records = list(SeqIO.parse(StringIO(fasta_text), "fasta")) #a list is made of the content

    # At least two sequences are required for comparison
    if len(records) < 2:
        raise ValueError("Please upload at least two FASTA sequences.")

    # Extract sequence identifiers
    ids = [record.id for record in records]

    # Extract DNA sequences and convert to uppercase
    sequences = [
        str(record.seq).upper()
        for record in records
    ]

    # Allowed nucleotide characters
    allowed = set("ACGTN-")

    # Identify sequences containing invalid characters
    invalid = [
        rid
        for rid, seq in zip(ids, sequences)
        if any(base not in allowed for base in seq) #checks every character (base) in seq. 
        #If even one character is not in allowed (such as a number, space, or unexpected letter), it returns True.
    ]

    if invalid:
        raise ValueError(
            "Unsupported characters in sequence(s): "
            + ", ".join(invalid)
        )

    # Check whether all sequences have equal lengths
    if len({len(seq) for seq in sequences}) != 1:
        raise ValueError(
            "Sequences have different lengths. "
            "Upload an aligned FASTA file."
        )

    return ids, sequences


# ------------------------------------------------------------
# PERSON B - FUNCTION 4: FASTA VARIANT ANALYSIS
# ------------------------------------------------------------

def analyse_fasta(sequences, reference):

    # Store identified variants
    rows = []

    # Compare every position in the reference sequence
    # enumerate() gives each reference base its Python position 
    for pos, ref_base in enumerate(reference):

        # Skip unknown reference bases and alignment gaps
        if ref_base not in "ACGT":
            continue

        # Collect valid nucleotides at the current position
        called = [
            seq[pos]
            for seq in sequences
            if seq[pos] in "ACGT"
        ]

        # Skip positions with no valid observations
        if not called:
            continue

        # Count nucleotide occurrences
        counts = Counter(called)

        # Calculate total valid observations
        denominator = len(called)

        # Identify nucleotides different from the reference
        for alt, count in sorted(counts.items()):

            if alt != ref_base:

                # Store variant information and frequency
                rows.append({
                    "Position (1-based)": pos + 1,
                    "Reference": ref_base,
                    "Alternate": alt,
                    "Alternate Count": count,
                    "Total Called Sequences": denominator,
                    "Frequency (%)": round(
                        count / denominator * 100, 2
                    )
                })

    # Convert results into a Pandas DataFrame
    return pd.DataFrame(rows)


# ============================================================
# PERSON A: STREAMLIT WEB APPLICATION
# ============================================================
# This section creates the actual website.
# It collects inputs, uploads files and displays results.
# ============================================================

# ------------------------------------------------------------
# PERSON A - SECTION 1: PAGE CONFIGURATION
# ------------------------------------------------------------

# Set browser-tab details and use the wide page layout.
st.set_page_config(
    page_title="Genomic Variant Frequency Calculator",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Genomic Variant Frequency Calculator")

st.caption(
    "Calculate variant frequencies from counts, CSV datasets, "
    "and aligned FASTA sequences."
)


# ------------------------------------------------------------
# PERSON A - SECTION 2: CREATE WEBSITE TABS
# ------------------------------------------------------------

# Create three tabs so each type of calculation has its own workspace.
tab1, tab2, tab3 = st.tabs([
    "Single Variant",
    "CSV Dataset",
    "FASTA Analysis"
])


# ------------------------------------------------------------
# PERSON A - SECTION 3: SINGLE VARIANT INTERFACE
# ------------------------------------------------------------

# Everything indented under 'with tab1' appears in the Single Variant tab.
with tab1:

    st.subheader("Calculate one variant")

    # text_input lets the user type a variant name/identifier.
    variant_id = st.text_input(
        "Variant ID",
        placeholder="e.g. rs123456"
    )

    variant_allele = st.text_input(
        "Variant allele",
        placeholder="e.g. A",
        max_chars=50
    )

    # number_input restricts the total to a positive whole-number count.
    total = st.number_input(
        "Total alleles observed",
        min_value=1,
        step=1,
        value=100
    )

    # This is the number of observations carrying the variant allele.
    variant_count = st.number_input(
        "Variant alleles observed",
        min_value=0,
        step=1,
        value=10
    )

    # The calculation runs only after the user clicks this button.
    if st.button(
        "Calculate frequency",
        type="primary",
        key="single_calc"
    ):

        try:
            result = calculate_single_frequency(
                variant_id,
                variant_allele,
                variant_count,
                total
            )

            st.success(
                f"Variant allele frequency: "
                f"**{result['Frequency (%)']:.2f}%**"
            )

            st.dataframe(
                pd.DataFrame([result]),
                use_container_width=True,
                hide_index=True
            )

        # If validation fails, show the reason on the page instead of
        # allowing the app to stop with an unhandled error.
        except Exception as exc:
            st.error(str(exc))
        
# ------------------------------------------------------------
# PERSON A - SECTION 4: CSV UPLOAD INTERFACE
# ------------------------------------------------------------

# Everything indented under 'with tab2' appears in the CSV Dataset tab.
with tab2:

    st.subheader("Calculate frequencies from a CSV file")

    st.markdown(
        "CSV columns required: "
        "`variant_id`, `variant_alleles`, `total_alleles`"
    )

    # file_uploader lets the user select a CSV from their device.
    csv_file = st.file_uploader(
        "Upload CSV",
        type=["csv"],
        key="csv_upload"
    )

    # A file uploader returns None until the user chooses a file.
    # Only try to process the CSV after a file has been uploaded.
    if csv_file is not None:

        try:
            df = process_csv(csv_file)

            st.success(
                f"Calculated frequencies for {len(df)} variant(s)."
            )

            st.dataframe(
                df,
                use_container_width=True,
                hide_index=True
            )

            # Convert the result table back into CSV text so the user
            # can download it. index=False avoids adding row numbers.
            st.download_button(
                "Download results as CSV",
                df.to_csv(index=False).encode("utf-8"),
                "variant_frequency_results.csv",
                "text/csv",
                key="csv_download"
            )
        except Exception as exc:
            st.error(f"Could not process CSV: {exc}")


# ------------------------------------------------------------
# PERSON A - SECTION 5: FASTA UPLOAD INTERFACE
# ------------------------------------------------------------

# Everything indented under 'with tab3' appears in the FASTA Analysis tab.
with tab3:

    st.subheader("Identify variants from aligned FASTA sequences")

    st.write(
        "Upload a multiple-sequence alignment in FASTA format. "
        "All sequences must already be aligned and have equal length."
    )

    st.info(
        "Each FASTA record is treated as one observed sequence/allele. "
        "The reported frequency is not automatically diploid clinical VAF."
    )

    # Accept common FASTA filename extensions in this upload control.
    fasta_file = st.file_uploader(
        "Upload FASTA file",
        type=["fasta", "fa", "fna", "fas"],
        key="fasta_upload"
    )

    # Do not parse anything until the user has uploaded a FASTA file.
    if fasta_file is not None:

        try:
            # Parse uploaded FASTA bytes using the corrected function.
            ids, sequences = parse_fasta(fasta_file.getvalue())

            ref_id = st.selectbox(
                "Select reference sequence",
                ids
            )

            # Find the selected reference's position in the IDs list,
            # then use the same position to get its DNA sequence.
            ref_index = ids.index(ref_id)
            reference = sequences[ref_index]

            # Run the comparison only when the user clicks Analyse.
            if st.button(
                "Analyse FASTA sequences",
                type="primary",
                key="fasta_analyse"
            ):

                result = analyse_fasta(
                    sequences,
                    reference
                )

                if not result.empty:

                    st.success(
                        f"Found {len(result)} alternate-allele result(s)."
                    )

                    st.dataframe(
                        result,
                        use_container_width=True,
                        hide_index=True
                    )

                    st.download_button(
                        "Download FASTA results",
                        result.to_csv(index=False).encode("utf-8"),
                        "fasta_variant_frequency_results.csv",
                        "text/csv",
                        key="fasta_download"
                    )

                else:
                    st.success(
                        "No nucleotide differences from the "
                        "selected reference were found."
                    )

            st.caption(
                f"Loaded {len(sequences)} sequences, "
                f"each {len(reference)} nucleotides long."
            )

        except Exception as exc:
            st.error(f"Could not process FASTA file: {exc}")


# ------------------------------------------------------------
# PERSON A - SECTION 6: INFORMATION DISPLAY
# ------------------------------------------------------------

# A divider separates the main tools from the explanatory information.
st.divider()

# This expandable panel explains the formulas and scientific limitations.
with st.expander("Calculation method and limitations"):

    st.markdown("""
    **Count-based frequency:**

    Variant alleles / Total alleles × 100.

    **FASTA-based alternate allele frequency:**

    Sequences carrying the alternate base /
    Sequences with a called base × 100.

    FASTA analysis uses Biopython `Bio.SeqIO` to read records.
    Sequences must be pre-aligned and equal in length.

    N and gaps (-) are excluded from the denominator.
    Only A, C, G and T are counted.

    Positions use 1-based numbering.""")

