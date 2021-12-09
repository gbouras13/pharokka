import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def run_phanotate(filepath):
    print("Beginning Phanotate")
    pack_dir = "prokka-env/bin/"
    out_dir = "output/"

    try:
        sp.call([os.path.join(pack_dir, "phanotate.py"), filepath, "-o", os.path.join(out_dir, "phanotate_out.fasta"), "-f", "fasta"], stderr=sp.DEVNULL, stdout=sp.DEVNULL) # silence the warnings (no trnaScan)
        sp.call([os.path.join(pack_dir, "phanotate.py"), filepath, "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"], stderr=sp.DEVNULL, stdout=sp.DEVNULL)
    except:
        sys.stderr.write("Error: phanotate not found")  
        return []

def translate_fastas():
    out_dir = "output/"
    with open(os.path.join(out_dir, "phanotate_aas.fasta"), 'w') as aa_fa:
        for dna_record in SeqIO.parse(os.path.join(out_dir, "phanotate_out.fasta"), 'fasta'): 
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_record.id)
            SeqIO.write(aa_record, aa_fa, 'fasta')
