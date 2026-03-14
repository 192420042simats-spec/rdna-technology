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

st.write("Upload plasmid FASTA sequence to analyze ORFs and possible MGE genes")

uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta","fa","txt"])


# ---------- CLEAN SEQUENCE ----------
def clean_sequence(seq):

    seq = seq.upper()

    seq = re.sub(r'[^ATGC]', '', seq)

    return seq


# ---------- ORF PREDICTION ----------
def predict_orfs(sequence):

    seq = Seq(sequence)

    proteins = []

    frames = [
        seq,
        seq[1:],
        seq[2:],
        seq.reverse_complement(),
        seq.reverse_complement()[1:],
        seq.reverse_complement()[2:]
    ]

    for frame in frames:

        trans = frame.translate()

        protein = ""

        for aa in str(trans):

            if aa == "*":

                if len(protein) > 50:
                    proteins.append(protein)

                protein = ""

            else:
                protein += aa

    return proteins


# ---------- FUNCTION PREDICTION ----------
def predict_function(protein):

    length = len(protein)

    if 300 <= length <= 450:
        return "Transposase (possible)"

    elif 250 <= length < 300:
        return "Integrase (possible)"

    elif 200 <= length < 250:
        return "Recombinase (possible)"

    elif 100 <= length < 200:
        return "Hypothetical protein"

    else:
        return "Unknown protein"


# ---------- PHYLOGENETIC TREE ----------
def build_phylogenetic_tree(proteins):

    max_len = max(len(p) for p in proteins)

    records = []

    for i, p in enumerate(proteins):

        padded = p.ljust(max_len, "-")

        records.append(
            SeqRecord(
                Seq(padded),
                id=f"ORF_{i+1}"
            )
        )

    alignment = MultipleSeqAlignment(records)

    calculator = DistanceCalculator("identity")

    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()

    tree = constructor.upgma(distance_matrix)

    return tree


# ---------- MAIN APP ----------
if uploaded_file:

    try:

        fasta_text = io.StringIO(uploaded_file.getvalue().decode("utf-8"))

        records = list(SeqIO.parse(fasta_text, "fasta"))

        raw_seq = str(records[0].seq)

        seq = clean_sequence(raw_seq)


        # STEP 1
        st.header("Step 1: Plasmid Information")

        st.write("Sequence Length:", len(seq))

        gc = gc_fraction(seq)

        st.write("GC Content:", round(gc*100,2), "%")


        # STEP 2
        st.header("Step 2: Automated ORF Prediction")

        proteins = predict_orfs(seq)

        st.write("Number of ORFs detected:", len(proteins))

        df_orf = pd.DataFrame({

            "ORF_ID":[f"ORF_{i+1}" for i in range(len(proteins))],
            "Protein_Length":[len(p) for p in proteins]

        })

        st.dataframe(df_orf)


        # STEP 3
        st.header("Step 3: MGE Detection + Function Prediction")

        results = []

        for i,p in enumerate(proteins):

            orf_id = f"ORF_{i+1}"

            function = predict_function(p)

            blast_link = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&QUERY={p}"

            results.append([
                orf_id,
                len(p),
                function,
                blast_link
            ])

        df_mge = pd.DataFrame(

            results,

            columns=[
                "ORF_ID",
                "Protein_Length",
                "Predicted_Function",
                "BLAST_Link"
            ]
        )

        st.dataframe(df_mge)


        # STEP 4
        st.header("Step 4: Phylogenetic Analysis")

        lengths = [len(p) for p in proteins]

        fig, ax = plt.subplots()

        ax.plot(range(len(lengths)), lengths, marker="o")

        ax.set_xlabel("ORF Index")

        ax.set_ylabel("Protein Length")

        ax.set_title("Protein Evolution Plot")

        st.pyplot(fig)


        st.subheader("Phylogenetic Tree")

        if len(proteins) >= 3:

            tree = build_phylogenetic_tree(proteins)

            fig_tree = plt.figure(figsize=(8,6))

            ax_tree = fig_tree.add_subplot(111)

            Phylo.draw(tree, axes=ax_tree, do_show=False)

            st.pyplot(fig_tree)

        else:

            st.write("Need at least 3 ORFs for tree.")


        # STEP 5
        st.header("Step 5: Genome Visualization")

        fig2, ax2 = plt.subplots()

        ax2.bar(range(len(lengths)), lengths)

        ax2.set_xlabel("Gene Index")

        ax2.set_ylabel("Gene Length")

        ax2.set_title("Plasmid Gene Distribution")

        st.pyplot(fig2)


        # DOWNLOAD
        st.header("Download Results")

        st.download_button(

            "Download Results",

            df_mge.to_csv(index=False),

            "mge_results.csv"

        )


    except Exception as e:

        st.error(f"Processing error: {e}")
