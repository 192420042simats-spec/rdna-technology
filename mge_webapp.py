import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt
import io

st.title("Mobile Genetic Element Analysis Web App")

st.write("Upload plasmid FASTA file to analyze ORFs and possible MGE genes")

uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta","fa","txt"])


# -------- ORF Detection Function --------
def find_orfs(sequence):

    orfs = []
    seq = Seq(sequence)

    frames = [seq,
              seq[1:],
              seq[2:],
              seq.reverse_complement(),
              seq.reverse_complement()[1:],
              seq.reverse_complement()[2:]]

    for frame in frames:

        trans = frame.translate()

        start = 0

        while True:

            start = trans.find("M", start)

            if start == -1:
                break

            stop = trans.find("*", start)

            if stop != -1:

                length = stop - start

                if length > 50:
                    protein = trans[start:stop]
                    orfs.append(str(protein))

            start += 1

    return orfs


# -------- MAIN APP --------
if uploaded_file is not None:

    try:

        # Convert uploaded binary file to text
        fasta_text = io.StringIO(uploaded_file.getvalue().decode("utf-8"))

        records = list(SeqIO.parse(fasta_text, "fasta"))

        if len(records) == 0:
            st.error("No FASTA sequence detected")
            st.stop()

        seq = str(records[0].seq)

        st.subheader("Basic Plasmid Information")

        length = len(seq)
        gc = gc_fraction(seq)

        st.write("Sequence Length:", length)
        st.write("GC Content:", round(gc * 100, 2), "%")

        # -------- ORF Detection --------
        st.subheader("ORF Detection")

        orfs = find_orfs(seq)

        st.write("Number of ORFs found:", len(orfs))

        if len(orfs) > 0:

            df = pd.DataFrame({
                "ORF_ID": [f"ORF_{i+1}" for i in range(len(orfs))],
                "Protein_Sequence": orfs
            })

            st.dataframe(df)

            st.download_button(
                label="Download ORF sequences",
                data=df.to_csv(index=False),
                file_name="orf_sequences.csv",
                mime="text/csv"
            )

        else:
            st.warning("No ORFs detected")

        # -------- MGE Keyword Suggestion --------
        st.subheader("Possible Mobile Genetic Element Proteins")

        st.write("Check these protein types using BLAST:")

        keywords = [
            "Transposase",
            "Integrase",
            "Recombinase",
            "Insertion sequence protein",
            "Phage protein"
        ]

        for k in keywords:
            st.write("-", k)

        # -------- Visualization --------
        st.subheader("Example Phylogenetic Tree Visualization")

        fig, ax = plt.subplots()

        ax.plot([1,2,3,4],[3,5,2,6], marker="o")

        ax.set_xlabel("Sequence")
        ax.set_ylabel("Distance")

        ax.set_title("Example Phylogenetic Tree")

        st.pyplot(fig)

    except Exception as e:
        st.error(f"Error processing file: {e}")
       
