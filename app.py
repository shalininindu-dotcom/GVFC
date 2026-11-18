import streamlit as st
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
from io import BytesIO

st.set_page_config(page_title="Genomic Variant Frequency Calculator", page_icon="🧬", layout="wide")

# =====================================================================
# LOOK AND FEEL  (this part is only styling - no science happens here)
# CSS = the language that controls colours and fonts on web pages.
# Theme: a blue "lab notebook" - graph paper, navy ink, blue highlights.
# =====================================================================
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Fraunces:ital,wght@0,600;0,800;1,400&family=IBM+Plex+Mono:wght@400;500&family=Inter:wght@400;500;600&display=swap');

/* 1. Pale-blue graph paper background (two sets of thin lines) */
.stApp {
    background-color: #EEF4FB;
    background-image:
        linear-gradient(rgba(47, 93, 158, 0.10) 1px, transparent 1px),
        linear-gradient(90deg, rgba(47, 93, 158, 0.10) 1px, transparent 1px);
    background-size: 32px 32px;
}
html, body, .stApp, p, label, span, li, button, input { font-family: 'Inter', sans-serif; }
h2, h3, p, label, li, span, .stMarkdown { color: #14284B; }

/* 2. The main content becomes a notebook page with a hard blue shadow */
.block-container {
    background: #FBFDFF; border: 1.5px solid #14284B; border-radius: 8px;
    padding: 2.5rem 3rem !important; margin-top: 2rem;
    box-shadow: 8px 8px 0 #2F5D9E;
}

/* 3. Title section */
.hero { border-bottom: 2px solid #14284B; padding-bottom: 20px; margin-bottom: 20px; }
.tag  { font-family: 'IBM Plex Mono', monospace; color: #2F5D9E; letter-spacing: 0.12em; font-size: 0.8rem; }
.hero h1 { font-family: 'Fraunces', serif; color: #14284B; font-size: 3.2rem; font-weight: 800; line-height: 1.05; margin: 8px 0; }
.hero h1 em { font-weight: 400; background: #BFE0FF; padding: 0 8px; }   /* blue highlighter pen */
.hero p  { color: #4C6485; font-size: 1.1rem; margin: 0; }

/* 4. "Fig. 1" strip of DNA letters */
.fig { background: #FFFFFF; border: 1.5px solid #14284B; border-radius: 6px; padding: 14px 20px; margin-bottom: 24px; }
.fig .cap { font-family: 'IBM Plex Mono', monospace; font-size: 0.75rem; color: #4C6485; }
.fig .seq { font-family: 'IBM Plex Mono', monospace; font-size: 1.3rem; letter-spacing: 0.35em; line-height: 1.7; white-space: nowrap; overflow-x: auto; color: #14284B; }
.fig .hit { background: #BFE0FF; outline: 1.5px solid #14284B; }

/* 5. Tabs */
button[data-baseweb="tab"] { color: #4C6485; font-size: 16px; font-weight: 600; }
button[data-baseweb="tab"][aria-selected="true"] { color: #14284B; }
div[data-baseweb="tab-highlight"] { background-color: #2F5D9E !important; }

/* 6. Buttons */
.stButton > button, .stDownloadButton > button {
    background: #2F5D9E; color: #FFFFFF; border: none;
    border-radius: 6px; font-weight: 600; padding: 0.5rem 1.6rem;
}
.stButton > button:hover, .stDownloadButton > button:hover { background: #244A82; color: #FFFFFF; }

/* 7. Input boxes */
div[data-baseweb="input"] { background: #F3F8FE !important; border: 1.5px solid #14284B !important; border-radius: 6px !important; }
input { background: transparent !important; color: #14284B !important; font-family: 'IBM Plex Mono', monospace !important; }
</style>
""", unsafe_allow_html=True)

# ---- Title section ----
st.markdown("""
<div class="hero">
  <h1>🧬 Genomic Variant <em>Frequency</em> Calculator</h1>
  <p>Calculate variant frequencies from counts, CSV datasets, and aligned FASTA sequences.</p>
</div>
<div class="fig">
  <div class="cap">Fig. 1 · aligned sequences, variant column highlighted</div>
  <div class="seq">
    ATGCCGTA<span class="hit">C</span>GGATCTA<br>
    ATGCCGTA<span class="hit">C</span>GGATCTA<br>
    ATGCCGTA<span class="hit">A</span>GGATCTA
  </div>
</div>
""", unsafe_allow_html=True)

# =====================================================================
# THE SCIENCE  (same logic as before)
# =====================================================================
tab1, tab2, tab3 = st.tabs(["Single Variant", "CSV Dataset", "FASTA Analysis"])

# ---------- TAB 1: one variant from counts ----------
with tab1:
    st.subheader("Calculate one variant")
    variant_id = st.text_input("Variant ID", placeholder="e.g. rs123456")
    variant_allele = st.text_input("Variant allele", placeholder="e.g. A", max_chars=50)
    total = st.number_input("Total alleles observed", min_value=1, step=1, value=100)
    variant_count = st.number_input("Variant alleles observed", min_value=0, step=1, value=10)

    if st.button("Calculate frequency", type="primary", key="single_calc"):
        allele = variant_allele.strip().upper()
        if not variant_id.strip():
            st.error("Please enter a variant ID.")
        elif variant_count > total:
            st.error("Variant alleles cannot be greater than total alleles.")
        elif allele and any(base not in "ACGTN" for base in allele):
            st.error("Allele should contain only A, C, G, T, or N.")
        else:
            if allele:
                Seq(allele)  # Biopython sequence object
            freq = variant_count / total * 100   # the core formula
            st.success(f"Variant allele frequency: **{freq:.2f}%**")
            st.dataframe(pd.DataFrame([{
                "Variant ID": variant_id.strip(),
                "Variant Allele": allele or "Not specified",
                "Variant Alleles": int(variant_count),
                "Total Alleles": int(total),
                "Frequency (%)": round(freq, 2)
            }]), use_container_width=True, hide_index=True)

# ---------- TAB 2: many variants from a CSV ----------
with tab2:
    st.subheader("Calculate frequencies from a CSV file")
    st.markdown("CSV columns required: `variant_id`, `variant_alleles`, `total_alleles`")
    csv_file = st.file_uploader("Upload CSV", type=["csv"], key="csv_upload")
    if csv_file is not None:
        try:
            df = pd.read_csv(csv_file)
            required = {"variant_id", "variant_alleles", "total_alleles"}
            missing = required - set(df.columns)
            if missing:
                st.error("Missing column(s): " + ", ".join(sorted(missing)))
            else:
                df["variant_alleles"] = pd.to_numeric(df["variant_alleles"], errors="coerce")
                df["total_alleles"] = pd.to_numeric(df["total_alleles"], errors="coerce")
                if df[["variant_alleles", "total_alleles"]].isna().any().any():
                    st.error("Allele counts must be numeric.")
                elif (df["total_alleles"] <= 0).any():
                    st.error("Total allele counts must be greater than zero.")
                elif (df["variant_alleles"] < 0).any():
                    st.error("Variant allele counts cannot be negative.")
                elif (df["variant_alleles"] > df["total_alleles"]).any():
                    st.error("Variant allele counts cannot exceed total allele counts.")
                else:
                    df["frequency_percent"] = (df["variant_alleles"] / df["total_alleles"] * 100).round(2)
                    st.success(f"Calculated frequencies for {len(df)} variant(s).")
                    st.dataframe(df, use_container_width=True, hide_index=True)
                    st.download_button("Download results as CSV", df.to_csv(index=False).encode("utf-8"),
                                       "variant_frequency_results.csv", "text/csv", key="csv_download")
        except Exception as exc:
            st.error(f"Could not read the CSV: {exc}")

# ---------- TAB 3: variants from an aligned FASTA file ----------
with tab3:
    st.subheader("Identify variants from aligned FASTA sequences")
    st.write("Upload a multiple-sequence alignment in FASTA format. All sequences must already be aligned and have equal length.")
    st.info("Each FASTA record is treated as one observed sequence/allele. The reported frequency is the proportion of valid sequences carrying an alternate base. It is not automatically diploid clinical VAF.")
    fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "fna", "fas"], key="fasta_upload")
    if fasta_file is not None:
        try:
            records = list(SeqIO.parse(BytesIO(fasta_file.getvalue()), "fasta"))
            if len(records) < 2:
                st.error("Please upload at least two FASTA sequences.")
            else:
                ids = [record.id for record in records]
                sequences = [str(record.seq).upper() for record in records]
                allowed = set("ACGTN-")
                invalid = [rid for rid, seq in zip(ids, sequences) if any(base not in allowed for base in seq)]
                if invalid:
                    st.error("Unsupported characters in sequence(s): " + ", ".join(invalid))
                elif len({len(seq) for seq in sequences}) != 1:
                    st.error("Sequences have different lengths. Upload an aligned FASTA file with equal-length sequences.")
                else:
                    ref_id = st.selectbox("Select reference sequence", ids)
                    ref_index = ids.index(ref_id)
                    reference = sequences[ref_index]
                    if st.button("Analyse FASTA sequences", type="primary", key="fasta_analyse"):
                        rows = []
                        for pos, ref_base in enumerate(reference):   # go through every position
                            if ref_base not in "ACGT":
                                continue
                            called = [seq[pos] for seq in sequences if seq[pos] in "ACGT"]
                            if not called:
                                continue
                            counts = Counter(called)                 # count how many of each base
                            denominator = len(called)
                            for alt, count in sorted(counts.items()):
                                if alt != ref_base:                  # only differences from the reference
                                    rows.append({
                                        "Position (1-based)": pos + 1,
                                        "Reference": ref_base,
                                        "Alternate": alt,
                                        "Alternate Count": count,
                                        "Total Called Sequences": denominator,
                                        "Frequency (%)": round(count / denominator * 100, 2)
                                    })
                        if rows:
                            result = pd.DataFrame(rows)
                            st.success(f"Found {len(result)} alternate-allele result(s).")
                            st.dataframe(result, use_container_width=True, hide_index=True)
                            st.download_button("Download FASTA results", result.to_csv(index=False).encode("utf-8"),
                                               "fasta_variant_frequency_results.csv", "text/csv", key="fasta_download")
                        else:
                            st.success("No nucleotide differences from the selected reference were found.")
                    st.caption(f"Loaded {len(records)} sequences, each {len(reference)} nucleotides long.")
        except Exception as exc:
            st.error(f"Could not read the FASTA file: {exc}")

st.divider()
with st.expander("Calculation method and limitations"):
    st.markdown("""
**Count-based frequency:** variant alleles / total alleles × 100.

**FASTA-based alternate allele frequency:** sequences carrying the alternate base / sequences with a called base at that position × 100.

FASTA analysis uses Biopython `Bio.SeqIO` to read records. Sequences must be pre-aligned and equal in length. `N` and gaps (`-`) are excluded from the denominator; only A, C, G and T are counted. Positions use 1-based numbering.

This educational tool compares aligned sequences. It does not map raw sequencing reads, call variants from raw reads, or determine clinical significance.
""")
st.caption("Built with Biopython • Streamlit • Pandas")
