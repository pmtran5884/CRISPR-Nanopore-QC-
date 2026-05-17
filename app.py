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
    st.markdown("**Note: Flanks must be present in your PCR amplicon!**")
    flank_5 = st.text_input("5' Flank (e.g., end of promoter)", value="CGAAACACCG").upper()
    flank_3 = st.text_input("3' Flank (e.g., start of scaffold)", value="GTTTAAGAGC").upper()
    min_q_score = st.slider("Minimum Average Q-Score", min_value=5, max_value=15, value=10)
    max_errors = st.number_input("Max Levenshtein Errors for Library Match", min_value=0, max_value=5, value=3)

# --- 2. Data Parsing & Processing ---
if fastq_file and key_file:
    if st.button("Run QC Analysis", type="primary"):
        
        df_key = pd.read_csv(key_file)
        library_guides = df_key['Guide_Seq'].str.upper().tolist()
        
        # --- Diagnostic Counters ---
        total_reads = 0
        passed_qc_reads = 0
        flanks_found = 0
        correct_length_barcodes = 0
        matched_to_library = 0
        
        extracted_barcodes = []
        all_extracted_lengths = []
        
        fastq_string = io.StringIO(fastq_file.getvalue().decode("utf-8"))
        
        with st.spinner("Loading FASTQ file into memory..."):
            all_reads = list(SeqIO.parse(fastq_string, "fastq"))
            total_reads = len(all_reads)
        
        progress_text = "Parsing reads and extracting barcodes..."
        my_bar = st.progress(0, text=progress_text)
        
        for i, record in enumerate(all_reads):
            if i % 1000 == 0 or i == total_reads - 1:
                progress_percent = int((i / total_reads) * 100)
                my_bar.progress(progress_percent, text=f"{progress_text} ({i}/{total_reads} reads)")
                
            # --- Filter 1: Quality Trimming ---
            quals = record.letter_annotations["phred_quality"]
            avg_q = sum(quals) / len(quals)
            if avg_q < min_q_score:
                continue
            passed_qc_reads += 1

            # --- Filter 2: Strand Handling & Flank Search ---
            seq_fwd = str(record.seq)
            seq_rev = str(record.seq.reverse_complement())
            
            b_fwd = ""
            b_rev = ""
            
            # Test Forward Strand (Added task="locations")
            m5_fwd = edlib.align(flank_5, seq_fwd, mode="HW", task="locations", k=2)
            m3_fwd = edlib.align(flank_3, seq_fwd, mode="HW", task="locations", k=2)
            if m5_fwd['editDistance'] != -1 and m3_fwd['editDistance'] != -1:
                start = m5_fwd['locations'][0][1] + 1
                end = m3_fwd['locations'][0][0]
                if start < end:
                    b_fwd = seq_fwd[start:end]

            # Test Reverse Strand (Added task="locations")
            m5_rev = edlib.align(flank_5, seq_rev, mode="HW", task="locations", k=2)
            m3_rev = edlib.align(flank_3, seq_rev, mode="HW", task="locations", k=2)
            if m5_rev['editDistance'] != -1 and m3_rev['editDistance'] != -1:
                start = m5_rev['locations'][0][1] + 1
                end = m3_rev['locations'][0][0]
                if start < end:
                    b_rev = seq_rev[start:end]

            # --- Filter 3: Valid Length Check ---
            barcode = None
            
            # Prioritize the strand that actually gives a ~20bp sequence
            if 17 <= len(b_fwd) <= 23:
                barcode = b_fwd
                flanks_found += 1
            elif 17 <= len(b_rev) <= 23:
                barcode = b_rev
                flanks_found += 1
            elif len(b_fwd) > 0 or len(b_rev) > 0:
                # Flanks were found, but the distance is garbage
                flanks_found += 1
                if len(b_fwd) > 0: all_extracted_lengths.append(len(b_fwd))
                elif len(b_rev) > 0: all_extracted_lengths.append(len(b_rev))
                
            if barcode:
                correct_length_barcodes += 1
                extracted_barcodes.append(barcode)

        my_bar.empty() 
        
        # --- Filter 4: Fuzzy Matching to Library ---
        guide_counts = {guide: 0 for guide in library_guides}
        
        with st.spinner(f'Matching {len(extracted_barcodes)} extracted barcodes to library...'):
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
                    matched_to_library += 1

        # --- Diagnostics Output UI ---
        st.success("Processing Complete!")
        st.header("Pipeline Diagnostics")
        
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("1. Total Reads", total_reads)
        col2.metric("2. Passed QC", passed_qc_reads)
        col3.metric("3. Flanks Found", flanks_found)
        col4.metric("4. Valid Length", correct_length_barcodes)
        col5.metric("5. Library Match", matched_to_library)
        
        st.divider()

        st.subheader("Diagnostic: Distance Between Flanks")
        if all_extracted_lengths:
            fig_len = px.histogram(x=all_extracted_lengths, nbins=100, 
                                   labels={'x': 'Distance between flanks (bp)', 'y': 'Number of Reads'}, 
                                   title="What is the actual length of the extracted barcodes?")
            fig_len.update_xaxes(range=[-50, 300]) 
            st.plotly_chart(fig_len)
        else:
            st.info("No invalid lengths to display.")
        
        st.divider()

        # --- Aggregation and Scoring ---
        df_results = pd.DataFrame(list(guide_counts.items()), columns=['Guide_Seq', 'Read_Count'])
        df_merged = pd.merge(df_results, df_key, on='Guide_Seq')
        
        gene_stats = df_merged.groupby('Gene').agg(
            Total_Reads=('Read_Count', 'sum'),
            Guides_Detected=('Read_Count', lambda x: (x > 0).sum()),
            Expected_Guides=('Guide_Seq', 'count')
        ).reset_index()
        
        gene_stats['Score'] = np.log2(1 + gene_stats['Total_Reads']) * (gene_stats['Guides_Detected'] / gene_stats['Expected_Guides'])
        gene_stats = gene_stats.sort_values(by='Score', ascending=False)

        # --- Output & Visualization ---
        st.header("Results Ranking")
        st.dataframe(gene_stats)
        
        csv = gene_stats.to_csv(index=False).encode('utf-8')
        st.download_button("Download Results CSV", csv, "crispr_qc_results.csv", "text/csv")
        
        st.header("Representation Quality")
        fig = px.scatter(gene_stats, x='Total_Reads', y='Guides_Detected', 
                         hover_data=['Gene'], title="Gene Representation vs. Read Depth")
        st.plotly_chart(fig)
