import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def run_phanotate(filepath):
    print("Beginning Phanotate")
    out_dir = "output/"
    try:
        sp.call([ "phanotate.py", filepath, "-o", os.path.join(out_dir, "phanotate_out.fasta"), "-f", "fasta"]) # , stderr=sp.DEVNULL, stdout=sp.DEVNULL silence the warnings (no trnaScan)
        sp.call(["phanotate.py", filepath, "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"])
    except:
        sys.stderr.write("Error: phanotate not found\n")  
        return 0

def tidy_phanotate_output():
    out_dir = "output/"
    phan_file = out_dir + "/phanotate_out.txt"
    col_list = ["start", "stop", "frame", "contig", "score"] 
    phan_df = pd.read_csv(phan_file, delimiter= '\t', index_col=False , names=col_list) 
    # get rid of the headers and reset the index
    phan_df = phan_df[phan_df['start'] != '#id:']
    phan_df = phan_df[phan_df['start'] != '#START'].reset_index(drop=True)
    phan_df["gene"] = phan_df['contig'] + " " + phan_df['start'] + "_" + phan_df['stop']
    phan_df.to_csv(out_dir + "cleaned_phanotate.tsv", sep="\t", index=False)
    return phan_df


def translate_fastas():
    print("Translating Nucleotide to Acids")
    phan_df = tidy_phanotate_output()
    out_dir = "output/"
    with open(os.path.join(out_dir, "phanotate_aas.fasta"), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, "phanotate_out.fasta"), 'fasta'): 
            dna_header = phan_df['contig'].iloc[i] 
            dna_description =   phan_df['start'].iloc[i] + "_" + phan_df['stop'].iloc[i]
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1

def run_trna_scan(filepath):
    print("Beginning tRNAscan-SE")
    out_dir = "output/"

    try:
        sp.call(["tRNAscan-SE", filepath, "-B", "-j",  os.path.join(out_dir, "trnascan_out.gff")])
    except:
        sys.stderr.write("Error: tRNAscan-SE not found\n")  
        return 0

    
              
def run_mmseqs():
    out_dir = "output/"
    phrog_db_dir = "databases/phrogs_mmseqs_db/"
    mmseqs_dir = "output/mmseqs/"
    amino_acid_fasta = "phanotate_aas.fasta"
    target_db_dir = "databases/target_dir"

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    sp.call(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")])
    # runs the seacr
    sp.call(["mmseqs", "search", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), "./tmp", "-s", "7"])
    sp.call(["mmseqs", "createtsv", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"mmseqs_results.tsv"), "--full-header"])
    # remove the target dir when finished 
    sp.call(["rm", "-r", target_db_dir])
