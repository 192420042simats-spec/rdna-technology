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

