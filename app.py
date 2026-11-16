import streamlit as st
import pandas as pd
from io import StringIO
from Bio.Seq import Seq

st.set_page_config(
    page_title="Genomic Variant Frequency Calculator",
    page_icon="🧬",
    layout="centered"
)

st.title("🧬 Genomic Variant Frequency Calculator")
st.caption("Calculate variant allele frequencies from genomic variant counts.")

st.markdown("""
This tool calculates **Variant Allele Frequency (VAF)** using:

**VAF (%) = Variant Alleles ÷ Total Alleles × 100**

Biopython is used to validate DNA alleles and perform basic sequence handling.
""")

tab1, tab2 = st.tabs(["Single Variant", "CSV Dataset"])

with tab1:
    st.subheader("Calculate one variant")

    variant_id = st.text_input("Variant ID", placeholder="e.g. rs123456")
    variant_allele = st.text_input(
        "Variant allele",
        placeholder="e.g. A",
        max_chars=50
    )
    total_alleles = st.number_input(
        "Total alleles observed",
        min_value=1,
        step=1,
        value=100
    )
    variant_alleles = st.number_input(
        "Variant alleles observed",
        min_value=0,
        step=1,
        value=10
    )

    if st.button("Calculate frequency", type="primary"):
        if not variant_id.strip():
            st.error("Please enter a variant ID.")
        elif variant_alleles > total_alleles:
            st.error("Variant alleles cannot be greater than total alleles.")
        else:
            allele = variant_allele.strip().upper()

            if allele:
                try:
                    Seq(allele)
                    if any(base not in "ACGTN" for base in allele):
                        st.error("Allele should contain only A, C, G, T, or N.")
                        st.stop()
                except Exception:
                    st.error("Invalid allele sequence.")
                    st.stop()

            frequency = (variant_alleles / total_alleles) * 100

            st.success(f"Variant allele frequency: **{frequency:.2f}%**")

            result = pd.DataFrame([{
                "Variant ID": variant_id.strip(),
                "Variant Allele": allele if allele else "Not specified",
                "Variant Alleles": variant_alleles,
                "Total Alleles": total_alleles,
                "Frequency (%)": round(frequency, 2)
            }])

            st.dataframe(result, use_container_width=True, hide_index=True)

with tab2:
    st.subheader("Calculate frequencies from a CSV file")

    st.markdown("""
Your CSV should contain these columns:

`variant_id`, `variant_alleles`, `total_alleles`

Example:

```text
variant_id,variant_alleles,total_alleles
rs123456,10,100
rs789012,25,200
rs345678,4,80
```
""")

    uploaded_file = st.file_uploader("Upload CSV", type=["csv"])

    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)

            required = {"variant_id", "variant_alleles", "total_alleles"}
            missing = required - set(df.columns)

            if missing:
                st.error(
                    "Missing required column(s): "
                    + ", ".join(sorted(missing))
                )
            else:
                df["variant_alleles"] = pd.to_numeric(
                    df["variant_alleles"], errors="coerce"
                )
                df["total_alleles"] = pd.to_numeric(
                    df["total_alleles"], errors="coerce"
                )

                if df[["variant_alleles", "total_alleles"]].isna().any().any():
                    st.error("Variant and total allele counts must be numeric.")
                elif (df["total_alleles"] <= 0).any():
                    st.error("Total allele counts must be greater than zero.")
                elif (df["variant_alleles"] < 0).any():
                    st.error("Variant allele counts cannot be negative.")
                elif (df["variant_alleles"] > df["total_alleles"]).any():
                    st.error(
                        "Variant allele count cannot be greater than total allele count."
                    )
                else:
                    df["frequency_percent"] = (
                        df["variant_alleles"] / df["total_alleles"] * 100
                    ).round(2)

                    st.success(f"Calculated frequencies for {len(df)} variant(s).")
                    st.dataframe(
                        df,
                        use_container_width=True,
                        hide_index=True
                    )

                    csv_output = df.to_csv(index=False).encode("utf-8")
                    st.download_button(
                        "Download results as CSV",
                        data=csv_output,
                        file_name="variant_frequency_results.csv",
                        mime="text/csv"
                    )
        except Exception as e:
            st.error(f"Could not read the file: {e}")

st.divider()

with st.expander("About the calculation"):
    st.write("""
    Variant Allele Frequency (VAF) represents the proportion of observed
    alleles that carry a particular variant.

    For example, if 15 of 100 observed alleles contain a variant:

    **VAF = 15 / 100 × 100 = 15%**

    This calculator reports the frequency from the counts supplied by the
    user. It does not determine whether a variant is harmful, pathogenic,
    or clinically significant.
    """)

st.caption("Biopython • Streamlit • Pandas")
