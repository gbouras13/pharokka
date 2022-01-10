import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def run_phanotate(filepath_in, out_dir):
    print("Beginning Phanotate")
    try:
        sp.call([ "phanotate.py", filepath_in, "-o", os.path.join(out_dir, "phanotate_out.fasta"), "-f", "fasta"]) # , stderr=sp.DEVNULL, stdout=sp.DEVNULL silence the warnings (no trnaScan)
        sp.call(["phanotate.py", filepath_in, "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"])
    except:
        sys.exit("Error: phanotate not found\n")  

def tidy_phanotate_output(out_dir):
    phan_file = os.path.join(out_dir, "phanotate_out.txt")
    col_list = ["start", "stop", "frame", "contig", "score"] 
    phan_df = pd.read_csv(phan_file, delimiter= '\t', index_col=False , names=col_list) 
    # get rid of the headers and reset the index
    phan_df = phan_df[phan_df['start'] != '#id:']
    phan_df = phan_df[phan_df['start'] != '#START'].reset_index(drop=True)
    phan_df["gene"] = ""
    # to match with hmms
    for index, row in phan_df.iterrows():
        row["gene"] = row['contig'] + str(index) + " " + row['start'] + "_" + row['stop']
    phan_df.to_csv(os.path.join(out_dir,"cleaned_phanotate.tsv"), sep="\t", index=False)
    return phan_df


def translate_fastas(out_dir):
    phan_df = tidy_phanotate_output(out_dir)
    with open(os.path.join(out_dir, "phanotate_aas.fasta"), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, "phanotate_out.fasta"), 'fasta'): 
            dna_header = phan_df['contig'].iloc[i] + str(i) 
            dna_description =   phan_df['start'].iloc[i] + "_" + phan_df['stop'].iloc[i]
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1

def run_trna_scan(filepath_in, out_dir):
    print("Beginning tRNAscan-SE")
    try:
        sp.call(["tRNAscan-SE", filepath_in, "-B", "-j",  os.path.join(out_dir, "trnascan_out.gff")])
    except:
        sys.stderr.write("Error: tRNAscan-SE not found\n")  
        return 0

    
def run_mmseqs(db_dir, out_dir):
    print("Running mmseqs")
    phrog_db_dir = os.path.join(db_dir, "phrogs_mmseqs_db/")
    mmseqs_dir = os.path.join(out_dir, "mmseqs/")
    amino_acid_fasta = "phanotate_aas.fasta"
    target_db_dir =  os.path.join(out_dir, "target_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    sp.call(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")])
    # runs the seacr
    sp.call(["mmseqs", "search", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), "./tmp", "-s", "8.5"])
    sp.call(["mmseqs", "createtsv", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"mmseqs_results.tsv"), "--full-header"])
    # remove the target dir when finished 
    sp.call(["rm", "-r", target_db_dir])

def run_hmmer(db_dir, out_dir):
    print("Running mmseqs")
    phrog_db_dir = os.path.join(db_dir, "phrogs_mmseqs_db/")
    mmseqs_dir = os.path.join(out_dir, "mmseqs/")
    amino_acid_fasta = "phanotate_aas.fasta"
    target_db_dir =  os.path.join(out_dir, "target_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    sp.call(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")])
    # runs the seacr
    sp.call(["mmseqs", "search", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), "./tmp", "-s", "8.5"])
    sp.call(["mmseqs", "createtsv", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"mmseqs_results.tsv"), "--full-header"])
    # remove the target dir when finished 
    sp.call(["rm", "-r", target_db_dir])

def run_hmmsuite(db_dir, out_dir):
    print("Running hmmsuite")
    hmmsuite_db_dir = os.path.join(db_dir, "phrogs_hhsuite_db/") 
    amino_acid_fasta = "phanotate_aas.fasta"
    target_db_dir =  os.path.join(out_dir, "hhsuite_target_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # custom index
    print(os.path.join(target_db_dir, 'hhsuite_tsv_file.ff{data,index}'))
    print(os.path.join(out_dir, amino_acid_fasta))

    # indexes the file 
    #can't pass curly brackets to subprocess so need to run using os
    cmd = 'ffindex_from_fasta -s ' + os.path.join(target_db_dir, 'hhsuite_tsv_file.ff{data,index}') + " " + os.path.join(out_dir, amino_acid_fasta)
    os.system(cmd)


    # runs
    sp.run(["hhblits_omp", '-i', os.path.join(target_db_dir, 'hhsuite_tsv_file'), '-d', os.path.join(hmmsuite_db_dir, "phrogs"), '-M', 'first', '-n', '1', '-o',os.path.join(target_db_dir, "results_your_seq_VS_phrogs"), '-blasttab', os.path.join(target_db_dir, "results_tsv_file")])
