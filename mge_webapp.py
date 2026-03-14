import streamlit as st
import io
import pandas as pd
import matplotlib.pyplot as plt
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

st.title("Mobile Genetic Element Analysis Web App")
st.write("Upload plasmid FASTA sequence to analyze ORFs and detect Mobile Genetic Elements")

# ---------------- Upload FASTA ----------------
uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta","fa","txt"])


# ---------------- Clean Sequence ----------------
def clean_sequence(seq):
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)
    return seq


# ---------------- ORF Prediction ----------------
def predict_orfs(sequence):

    start_codons = ["ATG"]
    stop_codons = ["TAA","TAG","TGA"]

    proteins = []
    orf_data = []

    seq_len = len(sequence)

    for frame in range(3):

        i = frame

        while i < seq_len - 3:

            codon = sequence[i:i+3]

            if codon in start_codons:

                for j in range(i, seq_len-3, 3):

                    stop = sequence[j:j+3]

                    if stop in stop_codons:

                        orf_seq = sequence[i:j+3]

                        protein = str(Seq(orf_seq).translate(to_stop=True))

                        if len(protein) > 30:

                            proteins.append(protein)

                            orf_data.append({
                                "Start": i,
                                "End": j+3,
                                "Length": len(orf_seq),
                                "Protein_length": len(protein)
                            })

                        i = j
                        break

            i += 3

    return proteins, pd.DataFrame(orf_data)


# ---------------- MGE Detection ----------------
def detect_mge(proteins):

    results = []

    for p in proteins:

        if len(p) > 100:

            results.append({
                "Protein_length": len(p),
                "Prediction": "Possible Mobile Genetic Element Protein"
            })

    return pd.DataFrame(results)


# ---------------- MAIN PROGRAM ----------------
if uploaded_file is not None:

    fasta_data = uploaded_file.read().decode("utf-8")

    record = SeqIO.read(io.StringIO(fasta_data), "fasta")

    sequence = clean_sequence(str(record.seq))


    # ---------------- Sequence Statistics ----------------
    st.subheader("Sequence Statistics")

    st.write("Sequence Length:", len(sequence))
    st.write("GC Content:", round(gc_fraction(sequence)*100,2), "%")


    # ---------------- ORF Detection ----------------
    st.subheader("Step 1: ORF Prediction")

    proteins, orf_df = predict_orfs(sequence)

    st.write("Number of ORFs detected:", len(orf_df))

    st.dataframe(orf_df)


    # ---------------- MGE Detection ----------------
    st.subheader("Step 2: Mobile Genetic Element Detection")

    mge_df = detect_mge(proteins)

    st.dataframe(mge_df)


    # ---------------- Phylogenetic Analysis ----------------
    st.subheader("Step 3: Phylogenetic Analysis")

    if len(proteins) > 1:

        selected = proteins[:5]

        max_len = max(len(p) for p in selected)

        records = []

        for i, p in enumerate(selected):

            padded = p.ljust(max_len, "-")

            records.append(
                SeqRecord(Seq(padded), id=f"Protein_{i}")
            )

        alignment = MultipleSeqAlignment(records)

        calculator = DistanceCalculator("identity")

        dm = calculator.get_distance(alignment)

        constructor = DistanceTreeConstructor()

        tree = constructor.nj(dm)

        fig = plt.figure(figsize=(6,4))

        Phylo.draw(tree, do_show=False)

        st.pyplot(fig)

    else:

        st.write("Not enough proteins for phylogenetic tree")


    # ---------------- Genome Visualization ----------------
    st.subheader("Step 4: Genome Visualization")

    if len(orf_df) > 0:

        fig2, ax = plt.subplots()

        ax.scatter(orf_df["Start"], orf_df["Length"])

        ax.set_xlabel("Genome Position")
        ax.set_ylabel("ORF Length")
        ax.set_title("Plasmid Gene Distribution")

        st.pyplot(fig2)


    # ---------------- Download Results ----------------
    st.subheader("Download Results")

    csv = orf_df.to_csv(index=False).encode()

    st.download_button(
        label="Download ORF Results CSV",
        data=csv,
        file_name="orf_results.csv",
        mime="text/csv"
        )
