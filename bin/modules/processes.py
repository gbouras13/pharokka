import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def run_phanotate(filepath_in, out_dir):
    print("Beginning Phanotate")
    add_delim_trim_fasta(filepath_in, out_dir)

    try:
        sp.run([ "phanotate.py", os.path.join(out_dir, "input_fasta_delim.fasta"), "-o", os.path.join(out_dir, "phanotate_out.fasta"), "-f", "fasta"], check=True) # , stderr=sp.DEVNULL, stdout=sp.DEVNULL silence the warnings (no trnaScan)
        sp.run(["phanotate.py", os.path.join(out_dir, "input_fasta_delim.fasta"), "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"], check=True)
    except:
        sys.exit("Error: phanotate not found\n")  

def add_delim_trim_fasta(filepath_in, out_dir):
    with open(os.path.join(out_dir, "input_fasta_delim.fasta"), 'w') as na_fa:
        for dna_record in SeqIO.parse(filepath_in, 'fasta'): 
            # trim to first 26 chars of fasta header if too long
            if len(dna_record.id) > 26:
                print("Trimming fasta headers to the first 26 characters.")
                dna_record.id = dna_record.id[0:25]
            # add delim
            dna_header = dna_record.id + str("delim") 
            dna_description = ""
            dna_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description)
            SeqIO.write(dna_record, na_fa, 'fasta')


def tidy_phanotate_output(out_dir):
    phan_file = os.path.join(out_dir, "phanotate_out.txt")
    col_list = ["start", "stop", "frame", "contig", "score"] 
    phan_df = pd.read_csv(phan_file, delimiter= '\t', index_col=False , names=col_list, skiprows=2 ) 
    # get rid of the headers and reset the index
    #phan_df = phan_df[phan_df['start'] != '#id:']
    # phan_df = phan_df[phan_df['start'] != '#START'].reset_index(drop=True)

    # to match with hmms
    #phan_df["gene"] = ""
    phan_df['gene'] = phan_df['contig'] + phan_df.index.astype(str) + " " + phan_df['start'].astype(str) + "_" + phan_df['stop'].astype(str)
    # old code
    # for index, row in phan_df.iterrows():
    #    # print(row['contig'] + str(index) + " " + str(row['start']) + "_" + str(row['stop']))
    #     row["gene"] = row['contig'] + str(index) + " " + str(row['start']) + "_" + str(row['stop'])
    phan_df.to_csv(os.path.join(out_dir,"cleaned_phanotate.tsv"), sep="\t", index=False)
    return phan_df


def translate_fastas(out_dir):
    phan_df = tidy_phanotate_output(out_dir)
    with open(os.path.join(out_dir, "phanotate_aas.fasta"), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, "phanotate_out.fasta"), 'fasta'): 
            dna_header = phan_df['contig'].iloc[i] + str(i) 
            dna_description =   phan_df['start'].iloc[i].astype(str) + "_" + phan_df['stop'].iloc[i].astype(str)
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1

def run_trna_scan(filepath_in, out_dir):
    print("Beginning tRNAscan-SE")
    try:
        sp.run(["tRNAscan-SE", filepath_in, "-B", "-Q", "-j",  os.path.join(out_dir, "trnascan_out.gff")], check=True)
    except:
        sys.stderr.write("Error: tRNAscan-SE not found\n")  
        return 0

    
def run_mmseqs(db_dir, out_dir, threads):
    print("Running mmseqs")
    phrog_db_dir = os.path.join(db_dir, "phrogs_mmseqs_db/")
    mmseqs_dir = os.path.join(out_dir, "mmseqs/")
    amino_acid_fasta = "phanotate_aas.fasta"
    target_db_dir =  os.path.join(out_dir, "target_dir/") 
    tmp_dir = os.path.join(out_dir, "tmp_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    sp.run(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")], check=True)
    # runs the seacrh
    sp.run(["mmseqs", "search", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), tmp_dir, "-s", "8.5",
    "--threads", threads], check=True)
    # creates the tsv
    sp.run(["mmseqs", "createtsv", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"mmseqs_results.tsv"), "--full-header", "--threads", threads], check=True)
    # remove the target dir when finished 
    sp.run(["rm", "-r", target_db_dir], check=True)


def run_hmmsuite(db_dir, out_dir, threads):
    print("Running hmmsuite")
    hmmsuite_db_dir = os.path.join(db_dir, "phrogs_hhsuite_db/")
    amino_acid_fasta = "phanotate_aas.fasta"
    target_db_dir =  os.path.join(out_dir, "hhsuite_target_dir/")
    tsv_prefix = os.path.join(target_db_dir, 'hhsuite_tsv_file.ff') 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)


    # indexes the file 
    sp.run(
        [
            'ffindex_from_fasta',
            '-s',
            ''.join((tsv_prefix, 'data')),
            ''.join((tsv_prefix, 'index')),
            os.path.join(out_dir, amino_acid_fasta),
        ],
        check=True
    )


    # runs
    sp.run(["hhblits_omp", '-i', os.path.join(target_db_dir, 'hhsuite_tsv_file'), '-d', os.path.join(hmmsuite_db_dir, "phrogs"), 
    '-M', 'first', '-n', '1', '-o',os.path.join(target_db_dir, "results_your_seq_VS_phrogs"), 
    '-blasttab', os.path.join(target_db_dir, "results_tsv_file"), "-cpu", threads], 
    check=True)
