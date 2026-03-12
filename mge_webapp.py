import streamlit as st
import io
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

st.title("Mobile Genetic Element Analysis Web App")

st.write("Upload plasmid DNA FASTA file")


uploaded_file = st.file_uploader(
    "Upload FASTA file",
    type=["fasta","fa","txt"]
)


# ---------------- ORF PREDICTION ----------------
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

        trans = frame.translate(to_stop=False)

        protein = ""

        for aa in str(trans):

            if aa == "*":
                if len(protein) > 50:
                    proteins.append(protein)
                protein = ""
            else:
                protein += aa

    return proteins


# ---------------- MGE DETECTION ----------------
def detect_mge(proteins):

    mge_keywords = [
        "transposase",
        "integrase",
        "recombinase",
        "insertion",
        "phage"
    ]

    results = []

    for i, p in enumerate(proteins):

        label = "Unknown"

        if len(p) > 300:
            label = "Possible Transposase"

        elif len(p) > 200:
            label = "Possible Integrase"

        results.append([f"ORF_{i+1}", len(p), label])

    return results


# ---------------- MAIN APP ----------------
if uploaded_file:

    try:

        fasta_text = io.StringIO(
            uploaded_file.getvalue().decode("utf-8")
        )

        records = list(SeqIO.parse(fasta_text, "fasta"))

        seq = str(records[0].seq).upper()

        st.header("Step 1: Basic Plasmid Information")

        st.write("Sequence Length:", len(seq))
        st.write("GC Content:", round(gc_fraction(seq)*100,2), "%")


        # ORF Prediction
        st.header("Step 2: Automated ORF Prediction")

        proteins = predict_orfs(seq)

        st.write("Number of ORFs detected:", len(proteins))

        df_orf = pd.DataFrame({
            "ORF_ID":[f"ORF_{i+1}" for i in range(len(proteins))],
            "Protein_Length":[len(p) for p in proteins]
        })

        st.dataframe(df_orf)


        # MGE Detection
        st.header("Step 3: MGE Detection")

        mge_results = detect_mge(proteins)

        df_mge = pd.DataFrame(
            mge_results,
            columns=["ORF_ID","Protein_Length","Predicted_Function"]
        )

        st.dataframe(df_mge)


        # Functional Annotation
        st.header("Step 4: Functional Annotation")

        st.write("Use BLAST to confirm function of detected ORFs.")


        # Phylogenetic Tree
        st.header("Step 5: Phylogenetic Analysis")

        fig, ax = plt.subplots()

        ax.plot(range(len(proteins)), [len(p) for p in proteins], marker="o")

        ax.set_xlabel("ORF")
        ax.set_ylabel("Protein Length")

        ax.set_title("Example Evolutionary Distance Plot")

        st.pyplot(fig)


        # Visualization
        st.header("Step 6: Genome Visualization")

        fig2, ax2 = plt.subplots()

        ax2.bar(range(len(proteins)), [len(p) for p in proteins])

        ax2.set_xlabel("Gene Index")
        ax2.set_ylabel("Gene Length")

        ax2.set_title("Gene Distribution on Plasmid")

        st.pyplot(fig2)

    except Exception as e:

        st.error(f"Processing error: {e}")
