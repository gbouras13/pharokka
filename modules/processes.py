import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def run_phanotate(filepath):
    print("Beginning Phanotate")
    out_dir = "output/"

    try:
        sp.call([ "phanotate.py", filepath, "-o", os.path.join(out_dir, "phanotate_out.fasta"), "-f", "fasta"]) # , stderr=sp.DEVNULL, stdout=sp.DEVNULL silence the warnings (no trnaScan)
        sp.call(["phanotate.py", filepath, "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"])
    except:
        sys.stderr.write("Error: phanotate not found\n")  
        return 0

def translate_fastas():
    print("Translating Nucleotide to Acids")
    out_dir = "output/"
    with open(os.path.join(out_dir, "phanotate_aas.fasta"), 'w') as aa_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "phanotate_out.fasta"), 'fasta'): 
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_record.id)
            SeqIO.write(aa_record, aa_fa, 'fasta')

def run_trna_scan(filepath):
    print("Beginning tRNAscan-SE")
    out_dir = "output/"

    try:
        sp.call(["tRNAscan-SE", filepath, "-B", "-o", os.path.join(out_dir, "trnascan_out.txt")])
    except:
        sys.stderr.write("Error: tRNAscan-SE not found\n")  
        return 0
              
def run_mmseqs():
    out_dir = "output/"
    db_dir = "databases/phrogs_mmseqs_db/"
    mmseqs_dir = "output/mmseqs/"
    amino_acid_fasta = "phanotate_aas.fasta"
    # creates db for input
    sp.call(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(out_dir, "target_seqs")])
    # runs the seacr
    sp.call(["mmseqs", "search", os.path.join(db_dir, "phrogs_profile_db"), os.path.join(out_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), "./tmp", "-s", "7"])
    sp.call(["mmseqs", "createtsv", os.path.join(db_dir, "phrogs_profile_db"), os.path.join(out_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"results.tsv"), "--full-header"])

