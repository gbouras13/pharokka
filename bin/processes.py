import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import logging
from BCBio import GFF
from Bio.Seq import Seq
from datetime import datetime

def write_to_log(s, logger):
           while True:
                output = s.readline().decode()
                if output:
                    logger.log(logging.INFO, output)
                else:
                    break


def run_phanotate(filepath_in, out_dir,logger):
    print("Running Phanotate.")
    logger.info("Running Phanotate.")
    try:
        phan_fast = sp.Popen(["phanotate.py", filepath_in, "-o", os.path.join(out_dir, "phanotate_out_tmp.fasta"), "-f", "fasta"], stderr=sp.PIPE, stdout=sp.DEVNULL) 
        phan_txt = sp.Popen(["phanotate.py", filepath_in, "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"], stderr=sp.PIPE, stdout=sp.DEVNULL)
        write_to_log(phan_fast.stderr, logger)
        write_to_log(phan_txt.stderr, logger)
    except:
        sys.exit("Error with Phanotate\n")  

def run_prodigal(filepath_in, out_dir,logger, meta, coding_table):
    print("Running Prodigal.")
    try:
        if meta == True:
            print("Prodigal Meta Mode Enabled.")
            logger.info("Prodigal Meta Mode Enabled.")
            prodigal = sp.Popen(["prodigal", "-i", filepath_in, "-d", os.path.join(out_dir, "prodigal_out_tmp.fasta"), "-f", "gff", "-o", os.path.join(out_dir, "prodigal_out.gff"), "-p", "meta", "-g", str(coding_table) ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        else:
            prodigal = sp.Popen(["prodigal", "-i", filepath_in, "-d", os.path.join(out_dir, "prodigal_out_tmp.fasta"), "-f", "gff", "-o", os.path.join(out_dir, "prodigal_out.gff"), "-g", str(coding_table) ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        write_to_log(prodigal.stdout, logger)
    except:
        sys.exit("Error with Prodigal\n")  


def tidy_phanotate_output(out_dir):
    phan_file = os.path.join(out_dir, "phanotate_out.txt")
    col_list = ["start", "stop", "frame", "contig", "score"] 
    phan_df = pd.read_csv(phan_file, delimiter= '\t', index_col=False , names=col_list, skiprows=2 ) 
    # get rid of the headers and reset the index
    phan_df = phan_df[phan_df['start'] != '#id:']
    phan_df = phan_df[phan_df['start'] != '#START'].reset_index(drop=True)
    phan_df['gene'] = phan_df['contig'] + phan_df.index.astype(str) + " " + phan_df['start'].astype(str) + "_" + phan_df['stop'].astype(str)
    phan_df.to_csv(os.path.join(out_dir,"cleaned_phanotate.tsv"), sep="\t", index=False)
    return phan_df

def tidy_prodigal_output(out_dir):
    prod_file = os.path.join(out_dir, "prodigal_out.gff")
    col_list = ["contig", "prod", "orf", "start", "stop","score", "frame", "phase", "description" ] 
    prod_df = pd.read_csv(prod_file, delimiter= '\t', index_col=False , names=col_list, skiprows=3 ) 
    
    # meta mode brings in some Nas so remove them
    prod_df = prod_df.dropna()
    prod_filt_df = prod_df[["start", "stop", "frame", "contig", "score"]]
    #convert staet stop to int
    prod_filt_df["start"] = prod_filt_df["start"].astype('int')
    prod_filt_df["stop"] = prod_filt_df["stop"].astype('int')
    # rearrange start and stop so that for negative strand, the stop is before start (like phanotate_out)
    cols = ["start","stop"]
    #indices where start is greater than stop
    ixs = prod_filt_df['frame'] == '-'
    # Where ixs is True, values are swapped
    prod_filt_df.loc[ixs,cols] = prod_filt_df.loc[ixs, cols].reindex(columns=cols[::-1]).values
    prod_filt_df['gene'] = prod_filt_df['contig'] + prod_filt_df.index.astype(str) + " " + prod_filt_df['start'].astype(str) + "_" + prod_filt_df['stop'].astype(str)
    prod_filt_df.to_csv(os.path.join(out_dir,"cleaned_prodigal.tsv"), sep="\t", index=False)
    return prod_filt_df


def translate_fastas(out_dir, gene_predictor):
    if gene_predictor == "phanotate":
        clean_df = tidy_phanotate_output(out_dir)
        fasta_input_tmp = "phanotate_out_tmp.fasta"
        fasta_output_aas_tmp = "phanotate_aas_tmp.fasta"
        fasta_output_aas_gd = "phanotate_aas.fasta"
        fasta_output_nts_gd = "phanotate_nts.fasta"
    if gene_predictor == "prodigal":
        clean_df = tidy_prodigal_output(out_dir)
        fasta_input_tmp = "prodigal_out_tmp.fasta"
        fasta_output_aas_tmp = "prodigal_aas_tmp.fasta"
        fasta_output_aas_gd = "prodigal_aas.fasta"
        fasta_output_nts_gd = "prodigal_nts.fasta"
    with open(os.path.join(out_dir, fasta_output_aas_tmp), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_tmp), 'fasta'): 
            dna_header = str(clean_df['contig'].iloc[i]) + str(i) 
            dna_description = str(clean_df['start'].iloc[i]) + "_" + str(clean_df['stop'].iloc[i])
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1
    with open(os.path.join(out_dir, fasta_output_aas_gd), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_tmp), 'fasta'): 
            dna_header = str(clean_df['contig'].iloc[i]) + "_" + str(i) 
            dna_description = str(clean_df['start'].iloc[i]) + "_" + str(clean_df['stop'].iloc[i])
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1
    with open(os.path.join(out_dir, fasta_output_nts_gd), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_tmp), 'fasta'): 
            dna_header = str(clean_df['contig'].iloc[i]) + "_" + str(i) 
            dna_description = str(clean_df['start'].iloc[i]) + "_" + str(clean_df['stop'].iloc[i])
            aa_record = SeqRecord(dna_record.seq, id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1
    

def run_trna_scan(filepath_in, out_dir, logger):
    print("Running tRNAscan-SE.")
    logger.info("Starting tRNA-scanSE")
    try:
        # needs stderr for trna scan
        trna = sp.Popen(["tRNAscan-SE", filepath_in, "-G", "-Q", "-j",  os.path.join(out_dir, "trnascan_out.gff")], stderr=sp.PIPE, stdout=sp.DEVNULL)
        write_to_log(trna.stderr, logger)
    except:
        sys.stderr.write("Error: tRNAscan-SE not found\n")  
        return 0

    
def run_mmseqs(db_dir, out_dir, threads, logger, gene_predictor, evalue):
    print("Running mmseqs2.")
    phrog_db_dir = os.path.join(db_dir, "phrogs_mmseqs_db/")
    mmseqs_dir = os.path.join(out_dir, "mmseqs/")
    amino_acid_fasta = gene_predictor + "_aas_tmp.fasta"
    target_db_dir =  os.path.join(out_dir, "target_dir/") 
    tmp_dir = os.path.join(out_dir, "tmp_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    mmseqs_createdb = sp.Popen(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")], stdout=sp.PIPE)
    write_to_log(mmseqs_createdb.stdout, logger)
    # runs the seacrh
    mmseqs_searc = sp.Popen(["mmseqs", "search", "-e", evalue ,os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), tmp_dir, "-s", "8.5",
    "--threads", threads], stdout=sp.PIPE)
    write_to_log(mmseqs_searc.stdout, logger)
    # creates the tsv
    mmseqs_createtsv = sp.Popen(["mmseqs", "createtsv", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"mmseqs_results.tsv"), "--full-header", "--threads", threads], stdout=sp.PIPE)
    write_to_log(mmseqs_createtsv.stdout, logger)
    # remove the target dir when finished 
    sp.run(["rm", "-r", target_db_dir], check=True)


def convert_gff_to_gbk(fasta_input, out_dir, prefix, logger):
    gff_file = os.path.join(out_dir, prefix + ".gff")
    gbk_file = os.path.join(out_dir, prefix + ".gbk")
    with open(gbk_file, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(fasta_input, "fasta"))
        for record in GFF.parse(gff_file, fasta_handler):
            for feature in record.features:
                # add translation only if CDS
                if feature.type == "CDS":
                    if feature.strand == 1:
                        feature.qualifiers.update({'translation': Seq.translate(record.seq[feature.location.start.position:feature.location.end.position], to_stop=True)})
                    else: # reverse strand -1 needs reverse compliment
                        feature.qualifiers.update({'translation': Seq.translate(record.seq[feature.location.start.position:feature.location.end.position].reverse_complement(), to_stop=True)})
                record.annotations["molecule_type"] = "DNA"
                record.annotations["date"] = datetime.today()
                record.annotations["topology"] = "linear"
                record.annotations["data_file_division"] = "VRL"
            SeqIO.write(record, gbk_handler, "genbank")

def run_minced(filepath_in, out_dir, prefix, logger):
    print("Running MinCED.")
    logger.info("Running MinCED.")
    try:
        # no phanotate stderr
        minced_fast = sp.Popen(["minced", filepath_in, os.path.join(out_dir, prefix + "_minced_spacers.txt") , os.path.join(out_dir, prefix + "_minced.gff")], stderr=sp.PIPE, stdout=sp.PIPE) 
        write_to_log(minced_fast.stderr, logger)
    except:
        sys.exit("Error with MinCED\n")  

def run_aragorn(filepath_in, out_dir, prefix, logger):
    print("Running Aragorn.")
    logger.info("Running Aragorn.")
    try:
        aragorn = sp.Popen(["aragorn", "-l", "-gcbact", "-w", "-o", os.path.join(out_dir, prefix + "_aragorn.txt"), "-m", filepath_in], stderr=sp.PIPE, stdout=sp.PIPE) 
        write_to_log(aragorn.stderr, logger)
    except:
        sys.exit("Error with Aragorn\n")  
