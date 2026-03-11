import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt

st.title("Mobile Genetic Element Analysis Web App")

st.write("Upload plasmid FASTA file to analyze ORFs and possible MGE genes")

uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta","fa","txt"])

def find_orfs(sequence):
    orfs = []
    seq = Seq(sequence)

    for frame in range(3):
        trans = seq[frame:].translate()

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


if uploaded_file:

    records = list(SeqIO.parse(uploaded_file, "fasta"))

    seq = str(records[0].seq)

    st.subheader("Basic Plasmid Information")

    length = len(seq)
    gc = gc_fraction(seq)

    st.write("Sequence Length:", length)
    st.write("GC Content:", round(gc*100,2), "%")

    st.subheader("ORF Detection")

    orfs = find_orfs(seq)

    st.write("Number of ORFs found:", len(orfs))

    if len(orfs) > 0:

        df = pd.DataFrame({
            "ORF_ID":[f"ORF_{i+1}" for i in range(len(orfs))],
            "Protein_Sequence":orfs
        })

        st.dataframe(df)

        st.download_button(
            label="Download ORF sequences",
            data=df.to_csv(index=False),
            file_name="orf_sequences.csv"
        )

    st.subheader("Possible Mobile Element Proteins")

    keywords = ["transposase","integrase","recombinase"]

    st.write("Check these proteins using BLAST for MGE identification")

    st.subheader("Phylogenetic Tree Placeholder")

    fig, ax = plt.subplots()
    ax.plot([1,2,3],[2,4,3])
    ax.set_title("Example Phylogenetic Tree Visualization")
    st.pyplot(fig)
