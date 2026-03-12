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

                            orf_data.append({
                                "ORF_ID": f"ORF_{len(proteins)}",
                                "Start": i,
                                "End": j+3,
                                "Length_nt": len(orf_seq),
                                "Protein_Length": len(protein),
                                "Protein_Sequence": protein
                            })

                        i = j
                        break

            i += 3

    return proteins, orf_data


# ---------------- MGE Detection ----------------
def predict_mge_type(protein):

    motifs = {
        "Transposase":["DDE","DDD"],
        "Integrase":["RHRY"],
        "Recombinase":["HNH"]
    }

    for mge, motif_list in motifs.items():
        for motif in motif_list:
            if motif in protein:
                return mge,"DNA Mobility"

    length = len(protein)

    if length > 300:
        return "Possible Transposase","Transposition"

    if length > 200:
        return "Possible Integrase","Integration"

    if length > 120:
        return "Insertion Sequence Protein","Mobility"

    return "Unknown Protein","Unknown"


# ---------------- Phylogenetic Tree ----------------
def build_tree(proteins):

    max_len = max(len(p) for p in proteins)

    records = []

    for i,p in enumerate(proteins):

        padded = p.ljust(max_len,"-")

        records.append(
            SeqRecord(Seq(padded), id=f"ORF_{i+1}")
        )

    alignment = MultipleSeqAlignment(records)

    calculator = DistanceCalculator("identity")

    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()

    tree = constructor.upgma(dm)

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


        # STEP 1
        st.header("Step 1: Plasmid Information")

        st.write("Sequence Length:",len(seq))

        gc = gc_fraction(seq)

        st.write("GC Content:",round(gc*100,2),"%")


        # STEP 2
        st.header("Step 2: Automated ORF Prediction")

        proteins, orf_data = predict_orfs(seq)

        st.write("Number of ORFs detected:",len(proteins))

        if len(proteins)==0:

            st.warning("No ORFs detected")

        df_orf = pd.DataFrame(orf_data)

        st.dataframe(df_orf)


        # STEP 3
        st.header("Step 3: MGE Identification")

        results = []

        for i,p in enumerate(proteins):

            mge,function = predict_mge_type(p)

            results.append([
                f"ORF_{i+1}",
                len(p),
                mge,
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


        # Highlight main ORF
        if len(df_mge)>0:

            main_orf = df_mge.loc[df_mge["Protein_Length"].idxmax()]

            st.subheader("Detected MGE ORF")

            st.success(f"MGE detected in {main_orf['ORF_ID']}")

            st.write("MGE Type:",main_orf["MGE_Type"])

            st.write("Function:",main_orf["Function"])

            st.code(main_orf["Protein_Sequence"])


        # STEP 4
        st.header("Step 4: Phylogenetic Analysis")

        lengths = [len(p) for p in proteins]

        fig,ax = plt.subplots()

        ax.plot(range(len(lengths)),lengths,marker="o")

        ax.set_xlabel("ORF Index")

        ax.set_ylabel("Protein Length")

        ax.set_title("Protein Evolution Plot")

        st.pyplot(fig)


        # Tree
        st.subheader("Phylogenetic Tree")

        if len(proteins) >= 3:

            try:

                tree = build_tree(proteins)

                fig_tree, ax_tree = plt.subplots(figsize=(8,6))

                Phylo.draw(tree, axes=ax_tree)

                st.pyplot(fig_tree)

            except:

                st.warning("Tree generation failed")

        else:

            st.warning("Need at least 3 ORFs to generate phylogenetic tree")


        # STEP 5
        st.header("Step 5: Genome Visualization")

        fig2,ax2 = plt.subplots()

        ax2.bar(range(len(lengths)),lengths)

        ax2.set_xlabel("Gene Index")

        ax2.set_ylabel("Gene Length")

        ax2.set_title("Plasmid Gene Distribution")

        st.pyplot(fig2)


        # STEP 6
        st.header("Download Results")

        csv = df_mge.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download MGE Results CSV",
            data=csv,
            file_name="mge_results.csv",
            mime="text/csv"
        )


    except Exception as e:

        st.error(f"Processing error: {e}")
