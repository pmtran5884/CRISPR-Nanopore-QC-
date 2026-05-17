import streamlit as st
import pandas as pd
import edlib
from Bio import SeqIO
from Bio.Seq import Seq
import plotly.express as px
import numpy as np
import io

# --- 1. User Inputs ---
st.title("Rapid CRISPR Screen QC (Nanopore)")

with st.sidebar:
    st.header("Upload Data")
    fastq_file = st.file_uploader("Upload Nanopore FASTQ (5-10MB)", type=['fastq', 'fq'])
    key_file = st.file_uploader("Upload sgRNA Key (CSV)", type=['csv'])
    
    st.header("Parameters")
    flank_5 = st.text_input("5' Flank (e.g., end of U6)", value="CGAAACACCG").upper()
    flank_3 = st.text_input("3' Flank (e.g., start of scaffold)", value="GTTTTAGAGC").upper()
    min_q_score = st.slider("Minimum Average Q-Score", min_value=5, max_value=15, value=10)
    max_errors = st.number_input("Max Levenshtein Errors for Library Match", min_value=0, max_value=5, value=3)

# --- 2. Data Parsing & Processing ---
if fastq_file and key_file:
    # Load Key into a dictionary for fast lookup
    df_key = pd.read_csv(key_file)
    library_guides = df_key['Guide_Seq'].str.upper().tolist()
    
    extracted_barcodes = []
    total_reads = 0
    passed_qc_reads = 0
    
    # Read FASTQ using io.StringIO for Streamlit UploadedFile compatibility
    fastq_string = io.StringIO(fastq_file.getvalue().decode("utf-8"))
    
    with st.spinner('Parsing and extracting reads...'):
        for record in SeqIO.parse(fastq_string, "fastq"):
            total_reads += 1
            
            # --- Quality Trimming ---
            quals = record.letter_annotations["phred_quality"]
            avg_q = sum(quals) / len(quals)
            if avg_q < min_q_score:
                continue
            
            passed_qc_reads += 1
            
            # --- Strand Handling ---
            # Test Forward Strand
            seq_fwd = str(record.seq)
            match_5 = edlib.align(flank_5, seq_fwd, mode="HW", k=2)
            match_3 = edlib.align(flank_3, seq_fwd, mode="HW", k=2)
            
            barcode = None
            if match_5['editDistance'] != -1 and match_3['editDistance'] != -1:
                start_idx = match_5['locations'][0][1] + 1
                end_idx = match_3['locations'][0][0]
                barcode = seq_fwd[start_idx:end_idx]
            else:
                # Test Reverse Complement Strand
                seq_rev = str(record.seq.reverse_complement())
                match_5_rev = edlib.align(flank_5, seq_rev, mode="HW", k=2)
                match_3_rev = edlib.align(flank_3, seq_rev, mode="HW", k=2)
                
                if match_5_rev['editDistance'] != -1 and match_3_rev['editDistance'] != -1:
                    start_idx = match_5_rev['locations'][0][1] + 1
                    end_idx = match_3_rev['locations'][0][0]
                    barcode = seq_rev[start_idx:end_idx]
            
            # Length Sanity Check (Allowing 17-23bp for a 20bp target due to indels)
            if barcode and 17 <= len(barcode) <= 23:
                extracted_barcodes.append(barcode)

    # --- 3. Fuzzy Matching to Library ---
    guide_counts = {guide: 0 for guide in library_guides}
    
    with st.spinner('Matching barcodes to library...'):
        for barcode in extracted_barcodes:
            best_match = None
            best_score = max_errors + 1 
            
            for guide in library_guides:
                dist = edlib.align(barcode, guide, mode="NW", task="distance")['editDistance']
                if dist < best_score:
                    best_score = dist
                    best_match = guide
                    
            if best_match:
                guide_counts[best_match] += 1

    # --- 4. Aggregation and Scoring ---
    st.success(f"Processing Complete! Passed QC: {passed_qc_reads}/{total_reads} reads.")
    
    df_results = pd.DataFrame(list(guide_counts.items()), columns=['Guide_Seq', 'Read_Count'])
    df_merged = pd.merge(df_results, df_key, on='Guide_Seq')
    
    # Group by Gene
    gene_stats = df_merged.groupby('Gene').agg(
        Total_Reads=('Read_Count', 'sum'),
        Guides_Detected=('Read_Count', lambda x: (x > 0).sum()),
        Expected_Guides=('Guide_Seq', 'count')
    ).reset_index()
    
    # Composite Score: Log2(reads) weighted by fraction of guides found
    gene_stats['Score'] = np.log2(1 + gene_stats['Total_Reads']) * (gene_stats['Guides_Detected'] / gene_stats['Expected_Guides'])
    gene_stats = gene_stats.sort_values(by='Score', ascending=False)

    # --- 5. Output & Visualization ---
    st.header("Results Ranking")
    st.dataframe(gene_stats)
    
    # Download Button
    csv = gene_stats.to_csv(index=False).encode('utf-8')
    st.download_button("Download Results CSV", csv, "crispr_qc_results.csv", "text/csv")
    
    st.header("Representation Quality")
    fig = px.scatter(gene_stats, x='Total_Reads', y='Guides_Detected', 
                     hover_data=['Gene'], title="Gene Representation vs. Read Depth")
    st.plotly_chart(fig)
