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
st.write("Upload a plasmid FASTA sequence to analyze ORFs and detect potential MGEs")

# ---------------- Upload FASTA ----------------
uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta","fa","txt"])


# ---------------- CLEAN SEQUENCE ----------------
def clean_sequence(seq):
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)
    return seq


# ---------------- ORF PREDICTION ----------------
def predict_orfs(sequence):

    start_codons = ["ATG"]
    stop_codons = ["TAA","TAG","TGA"]

    seq_len = len(sequence)
    proteins = []
    orf_info = []

    for frame in range(3):

        i = frame
        while i < seq_len-3:

            codon = sequence[i:i+3]

            if codon in start_codons:

                for j in range(i+3, seq_len-3, 3):

                    stop = sequence[j:j+3]

                    if stop in stop_codons:

                        orf_seq = sequence[i:j+3]

                        if len(orf_seq) > 150:

                            protein = str(Seq(orf_seq).translate(to_stop=True))

                            proteins.append(protein)

                            orf_info.append({
                                "Start": i,
                                "End": j+3,
                                "Length_nt": len(orf_seq),
                                "Protein_Length": len(protein),
                                "Protein_Sequence": protein
                            })

                        i = j
                        break
            i += 3

    return proteins, orf_info


# ---------------- MGE MOTIF DETECTION ----------------
def predict_mge_type(protein):

    motifs = {
        "Transposase": ["DDE","DDD","DEK"],
        "Integrase": ["RHRY","YREK"],
        "Recombinase": ["HNH","RHR"],
    }

    for mge, motif_list in motifs.items():
        for motif in motif_list:
            if motif in protein:
                return mge, "DNA mobility"

    length = len(protein)

    if length > 300:
        return "Possible Transposase", "Transposition"

    if length > 200:
        return "Possible Integrase", "Integration"

    if length > 120:
        return "Insertion Sequence Protein", "Mobility"

    return "Unknown Protein", "Unknown Function"


# ---------------- PHYLOGENETIC TREE ----------------
def build_phylogenetic_tree(proteins):

    max_len = max(len(p) for p in proteins)

    records = []

    for i,p in enumerate(proteins):

        padded = p.ljust(max_len,"-")

        records.append(
            SeqRecord(Seq(padded), id=f"ORF_{i+1}")
        )

    alignment = MultipleSeqAlignment(records)

    calculator = DistanceCalculator("identity")

    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()

    tree = constructor.upgma(distance_matrix)

    return tree


# ---------------- MAIN APP ----------------
if uploaded_file:

    try:

        fasta_text = io.StringIO(uploaded_file.getvalue().decode("utf-8"))

        records = list(SeqIO.parse(fasta_text,"fasta"))

        if len(records)==0:
            st.error("No sequence detected")
            st.stop()

        raw_seq = str(records[0].seq)

        seq = clean_sequence(raw_seq)

        # ---------------- STEP 1 ----------------
        st.header("Step 1: Plasmid Information")

        st.write("Sequence Length:",len(seq))

        gc = gc_fraction(seq)

        st.write("GC Content:",round(gc*100,2),"%")



        # ---------------- STEP 2 ----------------
        st.header("Step 2: ORF Prediction")

        proteins, orf_info = predict_orfs(seq)

        st.write("Number of ORFs detected:",len(proteins))

        if len(proteins)==0:

            st.warning("No ORFs detected")

            st.stop()

        df_orf = pd.DataFrame(orf_info)

        st.dataframe(df_orf)



        # ---------------- STEP 3 ----------------
        st.header("Step 3: MGE Identification")

        results = []

        for i,p in enumerate(proteins):

            mge_type,function = predict_mge_type(p)

            results.append([
                f"ORF_{i+1}",
                len(p),
                mge_type,
                function,
                p
            ])

        df_mge = pd.DataFrame(
            results,
            columns=[
                "ORF_ID",
                "Protein_Length",
                "MGE_Type",
                "Function",
                "Protein_Sequence"
            ]
        )

        st.dataframe(df_mge)



        # Highlight largest ORF
        main_orf = df_mge.loc[df_mge["Protein_Length"].idxmax()]

        st.subheader("Detected Major ORF")

        st.success(f"Main ORF: {main_orf['ORF_ID']}")

        st.write("Predicted Type:",main_orf["MGE_Type"])

        st.write("Function:",main_orf["Function"])

        st.code(main_orf["Protein_Sequence"])



        # ---------------- STEP 4 ----------------
        st.header("Step 4: Phylogenetic Analysis")

        lengths = [len(p) for p in proteins]

        fig,ax = plt.subplots()

        ax.plot(range(len(lengths)),lengths,marker="o")

        ax.set_xlabel("ORF Index")

        ax.set_ylabel("Protein Length")

        ax.set_title("Protein Length Evolution")

        st.pyplot(fig)



        # Tree
        st.subheader("Phylogenetic Tree")

        if len(proteins) >= 3:

            tree = build_phylogenetic_tree(proteins)

            fig_tree = plt.figure(figsize=(8,6))

            ax_tree = fig_tree.add_subplot(111)

            Phylo.draw(tree, axes=ax_tree, do_show=False)

            st.pyplot(fig_tree)

        else:

            st.info("Need at least 3 ORFs for tree")



        # ---------------- STEP 5 ----------------
        st.header("Step 5: Genome Visualization")

        fig2,ax2 = plt.subplots()

        ax2.bar(range(len(lengths)),lengths)

        ax2.set_xlabel("Gene Index")

        ax2.set_ylabel("Gene Length")

        ax2.set_title("Gene Distribution")

        st.pyplot(fig2)



        # ---------------- STEP 6 ----------------
        st.header("Download Results")

        st.download_button(
            "Download CSV",
            df_mge.to_csv(index=False),
            "MGE_results.csv"
        )



    except Exception as e:

        st.error(f"Processing Error: {e}")
