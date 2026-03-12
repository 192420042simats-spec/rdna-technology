import streamlit as st
import io
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

st.title("Mobile Genetic Element Analysis Web App")
st.write("Upload plasmid DNA FASTA sequence for ORF and MGE analysis")

uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta","fa","txt"])

# ---------------- VALIDATE DNA ----------------
def is_dna(seq):
    valid = set("ATGC")
    seq = seq.upper().replace("\n","").replace(" ","")
    return set(seq).issubset(valid)

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
    results = []
    for i, p in enumerate(proteins):
        length = len(p)
        if 300 <= length <= 450:
            label = "Possible Transposase"
        elif 250 <= length < 300:
            label = "Possible Integrase"
        elif 200 <= length < 250:
            label = "Possible Recombinase"
        elif length > 450:
            label = "Large Mobile Element Protein"
        else:
            label = "Unknown"
        results.append([f"ORF_{i+1}", length, label, p])
    return results

# ---------------- MAIN APP ----------------
if uploaded_file:
    try:
        fasta_text = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        records = list(SeqIO.parse(fasta_text, "fasta"))
        if len(records) == 0:
            st.error("No FASTA sequence detected")
            st.stop()
        seq = str(records[0].seq).upper().replace("\n","").replace(" ","")
        if not is_dna(seq):
            st.error("Uploaded file is not DNA FASTA. Please upload plasmid DNA sequence containing only A,T,G,C.")
            st.stop()

        # Step 1: Basic Info
        st.header("Step 1: Basic Plasmid Information")
        st.write("Sequence Length:", len(seq))
        st.write("GC Content:", round(gc_fraction(seq)*100,2), "%")

        # Step 2: ORF Prediction
        st.header("Step 2: Automated ORF Prediction")
        proteins = predict_orfs(seq)
        st.write("Number of ORFs detected:", len(proteins))
        df_orf = pd.DataFrame({
            "ORF_ID":[f"ORF_{i+1}" for i in range(len(proteins))],
            "Protein_Length":[len(p) for p in proteins]
        })
        st.dataframe(df_orf)

        # Step 3: MGE Detection
        st.header("Step 3: MGE Detection")
        mge_results = detect_mge(proteins)
        df_mge = pd.DataFrame(mge_results, columns=["ORF_ID","Protein_Length","Predicted_Function","Protein_Sequence"])
        st.dataframe(df_mge)

        # Step 4: Functional Annotation
        st.header("Step 4: Functional Annotation")
        st.write("Check these ORFs on NCBI BLAST to confirm function:")
        for i, row in df_mge.iterrows():
            protein = row["Protein_Sequence"]
            blast_link = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&QUERY={protein}"
            st.write(row["ORF_ID"], "→", blast_link)

        # Step 5: Phylogenetic Analysis
        st.header("Step 5: Phylogenetic Analysis")
        lengths = [len(p) for p in proteins]
        fig, ax = plt.subplots()
        ax.plot(range(len(lengths)), lengths, marker="o")
        ax.set_xlabel("ORF Index")
        ax.set_ylabel("Protein Length")
        ax.set_title("Protein Length Distance Plot")
        st.pyplot(fig)

        # Step 6: Genome Visualization
        st.header("Step 6: Genome Visualization")
        fig2, ax2 = plt.subplots()
        ax2.bar(range(len(lengths)), lengths)
        ax2.set_xlabel("Gene Index")
        ax2.set_ylabel("Gene Length")
        ax2.set_title("Plasmid Gene Distribution")
        st.pyplot(fig2)

        # Download results
        st.header("Download Results")
        st.download_button("Download ORF and MGE table", df_mge.to_csv(index=False), "mge_results.csv")

    except Exception as e:
        st.error(f"Processing error: {e}")
