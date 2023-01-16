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
import pyrodigal


def write_to_log(s, logger):
           while True:
                output = s.readline().decode()
                if output:
                    logger.log(logging.INFO, output)
                else:
                    break

##### phanotate meta mode ########

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    https://biopython.org/wiki/Split_large_file

    :param iterator: iterator for enumerating over
    :param batch_size: number of fasta records in each file
    :return:
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []

def split_input_fasta(filepath_in, out_dir):
    """Splits the input fasta into separate single fasta files for multithreading with phanotate
    https://biopython.org/wiki/Split_large_file

    :param filepath_in: input multifasta file
    :param out_dir: output director 
    :return: num_fastas: int giving the number of fasta records in the multifasta
    """
    # iterate and count fastas
    record_iter = SeqIO.parse(open(filepath_in), "fasta")
    num_fastas = len([1 for line in open(filepath_in) if line.startswith(">")])

    # each fasta gets its own file so batch size of 1
    batch_size = 1

    input_tmp_dir = os.path.join(out_dir, "input_split_tmp")

    for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
        filename = "input_subprocess%i.fasta" % (i + 1)
        with open(os.path.join(input_tmp_dir, filename), "w") as handle:
            SeqIO.write(batch, handle, "fasta")
    return num_fastas



def run_phanotate_fasta_meta(filepath_in, out_dir, threads, num_fastas):
    """
    Runs phanotate to output fastas
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param threads: threads
    :param num_fastas: number of fastas in input multifasta
    :return:
    """

    phanotate_tmp_dir = os.path.join(out_dir, "input_split_tmp")
    commands = []

    for i in range(1,num_fastas+1):
        in_file = "input_subprocess" + str(i) +".fasta" 
        out_file = "phanotate_out_tmp" + str(i) +".fasta"
        filepath_in = os.path.join(phanotate_tmp_dir,in_file)
        cmd = "phanotate.py " + filepath_in + " -o " + os.path.join(phanotate_tmp_dir, out_file) + " -f fasta"
        commands.append(cmd)

    n = int(threads) #the number of parallel processes you want

    for j in range(max(int(len(commands)/n)+1, 1)):
        procs = [sp.Popen(i, shell=True) for i in commands[j*n: min((j+1)*n, len(commands))] ]
        for p in procs:
            p.wait()


def run_phanotate_txt_meta(filepath_in, out_dir, threads, num_fastas):
    """
    Runs phanotate to output text file
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param threads: threads
    :param num_fastas: number of fastas in input multifasta
    :return:
    """

    phanotate_tmp_dir = os.path.join(out_dir, "input_split_tmp")

    commands = []

    for i in range(1,num_fastas+1):
        in_file = "input_subprocess" + str(i) +".fasta" 
        out_file = "phanotate_out_tmp" + str(i) +".txt"
        filepath_in = os.path.join(phanotate_tmp_dir,in_file)
        cmd = "phanotate.py " + filepath_in + " -o " + os.path.join(phanotate_tmp_dir, out_file) + " -f tabular"
        commands.append(cmd)

    n = int(threads) #the number of parallel processes you want
    for j in range(max(int(len(commands)/n)+1, 1)):
        procs = [sp.Popen(i, shell=True) for i in commands[j*n: min((j+1)*n, len(commands))] ]
        for p in procs:
            p.wait()

def concat_phanotate_meta(out_dir, num_fastas):
    """
    Concatenates phanotate output for downstream analysis
    :param out_dir: output directory
    :param threads: threads
    :return:
    """

    phanotate_tmp_dir = os.path.join(out_dir, "input_split_tmp")

    tsvs = []
    for i in range(1,int(num_fastas)+1):
        out_tsv = "phanotate_out_tmp" + str(i) +".txt"
        tsvs.append(os.path.join(phanotate_tmp_dir, out_tsv))

    with open(os.path.join(out_dir, "phanotate_out.txt"), 'w') as outfile:
        for fname in tsvs:
            with open(fname) as infile:
                outfile.write(infile.read())

    fastas = []
    for i in range(1,int(num_fastas)+1):
        out_fasta = "phanotate_out_tmp" + str(i) +".fasta"
        fastas.append(os.path.join(phanotate_tmp_dir, out_fasta))

    with open(os.path.join(out_dir, "phanotate_out_tmp.fasta"), 'w') as outfile:
        for fname in fastas:
            with open(fname) as infile:
                outfile.write(infile.read())


def run_trnascan_meta(filepath_in, out_dir, threads, num_fastas):
    """
    Runs trnascan to output gffs one contig per thread
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param threads: threads
    :param num_fastas: number of fastas in input multifasta
    :return:
    """

    input_tmp_dir = os.path.join(out_dir, "input_split_tmp")
    commands = []

    for i in range(1,num_fastas+1):
        in_file = "input_subprocess" + str(i) +".fasta" 
        out_file = "trnascan_tmp" + str(i) +".gff"
        filepath_in = os.path.join(input_tmp_dir,in_file)
        filepath_out = os.path.join(input_tmp_dir,out_file)
        cmd = "tRNAscan-SE " + filepath_in + " --thread 1 -G -Q -j " + filepath_out
        commands.append(cmd)

    n = int(threads) #the number of parallel processes you want

    for j in range(max(int(len(commands)/n)+1, 1)):
        procs = [sp.Popen(i, shell=True , stderr=sp.PIPE, stdout=sp.DEVNULL) for i in commands[j*n: min((j+1)*n, len(commands))] ]
        for p in procs:
            p.wait()

def concat_trnascan_meta(out_dir, num_fastas):
    """
    Concatenates trnascan output for downstream analysis
    :param out_dir: output directory
    :param threads: threads
    :return:
    """

    input_tmp_dir = os.path.join(out_dir, "input_split_tmp")

    gffs = []
    for i in range(1,int(num_fastas)+1):
        out_gff = "trnascan_tmp" + str(i) +".gff"
        gffs.append(os.path.join(input_tmp_dir, out_gff))

    with open(os.path.join(out_dir, "trnascan_out.gff"), 'w') as outfile:
        for fname in gffs:
            with open(fname) as infile:
                outfile.write(infile.read())

 



##### single contig mode ######

def run_phanotate(filepath_in, out_dir,logger):
    """
    Runs phanotate
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :return:
    """
    try:
        phan_fast = sp.Popen(["phanotate.py", filepath_in, "-o", os.path.join(out_dir, "phanotate_out_tmp.fasta"), "-f", "fasta"], stderr=sp.PIPE, stdout=sp.DEVNULL) 
        phan_txt = sp.Popen(["phanotate.py", filepath_in, "-o", os.path.join(out_dir, "phanotate_out.txt"), "-f", "tabular"], stderr=sp.PIPE, stdout=sp.DEVNULL)
        write_to_log(phan_fast.stderr, logger)
        write_to_log(phan_txt.stderr, logger)
    except:
        sys.exit("Error with Phanotate\n")  


def run_prodigal(filepath_in, out_dir,logger, meta, coding_table):
    """
    Gets CDS using prodigal
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :param meta Boolean - metagenomic mode flag 
    :param coding_table coding table for prodigal (default 11)
    :return:
    """
    print("Running Prodigal")
    try:
        if meta == True:
            print("Prodigal Meta Mode Enabled")
            logger.info("Prodigal Meta Mode Enabled")
            prodigal = sp.Popen(["prodigal", "-i", filepath_in, "-d", os.path.join(out_dir, "prodigal_out_tmp.fasta"), "-f", "gff", "-o", os.path.join(out_dir, "prodigal_out.gff"), "-p", "meta", "-g", str(coding_table) ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        else:
            prodigal = sp.Popen(["prodigal", "-i", filepath_in, "-d", os.path.join(out_dir, "prodigal_out_tmp.fasta"), "-f", "gff", "-o", os.path.join(out_dir, "prodigal_out.gff"), "-g", str(coding_table) ], stdout=sp.PIPE, stderr=sp.DEVNULL) 
        write_to_log(prodigal.stdout, logger)
    except:
        sys.exit("Error with Prodigal\n")  


def run_pyrodigal(filepath_in, out_dir,logger, meta, coding_table):
    """
    Gets CDS using pyrodigal
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :param meta Boolean - metagenomic mode flag 
    :param coding_table coding table for prodigal (default 11)
    :return:
    """

    prodigal_metamode = False
    if meta == True:
        prodigal_metamode = True 
        print("Prodigal Meta Mode Enabled")
        logger.info("Prodigal Meta Mode Enabled")

    # for training if you want different coding table
    seqs = [bytes(record.seq) for record in SeqIO.parse(filepath_in, 'fasta' )]
    record = SeqIO.parse(filepath_in, 'fasta' )
    orf_finder = pyrodigal.OrfFinder(meta=prodigal_metamode )

    # coding table possible if false
    if prodigal_metamode == False:
        trainings_info = orf_finder.train(*seqs, translation_table=int(coding_table))
        orf_finder = pyrodigal.OrfFinder(trainings_info, meta=prodigal_metamode )

    with open(os.path.join(out_dir, "prodigal_out.gff"), "w") as dst:
        for i, record in enumerate(SeqIO.parse(filepath_in, 'fasta' )):
            genes = orf_finder.find_genes(str(record.seq))
            genes.write_gff(dst, sequence_id=record.id)
    
    with open(os.path.join(out_dir, "prodigal_out_tmp.fasta"), "w") as dst:
        for i, record in enumerate(SeqIO.parse(filepath_in, 'fasta' )):
            genes = orf_finder.find_genes(str(record.seq))
            genes.write_genes(dst, sequence_id=record.id)

    





def tidy_phanotate_output(out_dir):
    """
    Tidies phanotate output
    :param out_dir: output directory
    :return: phan_df pandas dataframe
    """
    phan_file = os.path.join(out_dir, "phanotate_out.txt")
    col_list = ["start", "stop", "frame", "contig", "score"] 
    phan_df = pd.read_csv(phan_file, delimiter= '\t', index_col=False , names=col_list, skiprows=2 ) 
    # get rid of the headers and reset the index
    phan_df = phan_df[phan_df['start'] != '#id:']
    phan_df = phan_df[phan_df['start'] != '#START'].reset_index(drop=True)
    phan_df['gene'] = phan_df['contig'].astype(str) + phan_df.index.astype(str) + " " + phan_df['start'].astype(str) + "_" + phan_df['stop'].astype(str)
    phan_df.to_csv(os.path.join(out_dir,"cleaned_phanotate.tsv"), sep="\t", index=False)
    return phan_df

def tidy_prodigal_output(out_dir):
    """
    Tidies prodigal output
    :param out_dir: output directory
    :return: prod_filt_df pandas dataframe
    """
    prod_file = os.path.join(out_dir, "prodigal_out.gff")
    col_list = ["contig", "prod", "orf", "start", "stop","score", "frame", "phase", "description" ] 
    prod_df = pd.read_csv(prod_file, delimiter= '\t', index_col=False , names=col_list, skiprows=3 ) 
    
    # meta mode brings in some Nas so remove them
    # need to reset index!!!! and drop, or else will cause rubbish results for metagenomics
    prod_df = prod_df.dropna().reset_index(drop=True)


    prod_filt_df = prod_df[["start", "stop", "frame", "contig", "score"]]

    #convert start stop to int
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
    """
    Translates input CDSs to amino acids. For now will use 11 translation table. Will get around to alternative coding later
    :param out_dir: output directory
    :param gene_predictor: phanotate or prodigal
    :return: 
    """
    if gene_predictor == "phanotate":
        clean_df = tidy_phanotate_output(out_dir)
    if gene_predictor == "prodigal":
        clean_df = tidy_prodigal_output(out_dir)
    
    fasta_input_tmp = gene_predictor + "_out_tmp.fasta"
    fasta_output_aas_tmp = gene_predictor + "_aas_tmp.fasta"

    # translate for temporary AA output
    with open(os.path.join(out_dir, fasta_output_aas_tmp), 'w') as aa_fa:
        i = 0 
        for dna_record in SeqIO.parse(os.path.join(out_dir, fasta_input_tmp), 'fasta'): 
            dna_header = str(clean_df['contig'].iloc[i]) + str(i) 
            dna_description = str(clean_df['start'].iloc[i]) + "_" + str(clean_df['stop'].iloc[i])
            aa_record = SeqRecord(dna_record.seq.translate(to_stop=True), id=dna_header, description = dna_description )
            SeqIO.write(aa_record, aa_fa, 'fasta')
            i += 1


def run_trna_scan(filepath_in,threads, out_dir, logger):
    """
    Runs trna scan
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :return:
    """
    try:
        # needs stderr for trna scan
        trna = sp.Popen(["tRNAscan-SE", filepath_in, "--thread",threads, "-G", "-Q", "-j",  os.path.join(out_dir, "trnascan_out.gff")], stderr=sp.PIPE, stdout=sp.DEVNULL)
        write_to_log(trna.stderr, logger)
    except:
        sys.stderr.write("Error: tRNAscan-SE not found\n")  
        return 0

    
def run_mmseqs(db_dir, out_dir, threads, logger, gene_predictor, evalue):
    """
    Runs mmseqs2 on phrogs 
    :param db_dir: database path
    :param out_dir: output directory
    :param logger: logger
    :params threads: threads
    :param gene_predictor: phanotate or prodigal
    :param evalue: evalue for mmseqs2
    :return:
    """
    print("Running MMseqs2 on PHROGs Database.")
    logger.info("Running MMseqs2 on PHROGs Database.")

    # declare directories - phrog_db_dir is now the db_dir
    phrog_db_dir = db_dir
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
    # runs the mmseqs seacrh
    mmseqs_search = sp.Popen(["mmseqs", "search", "-e", evalue ,os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), tmp_dir, "-s", "8.5",
    "--threads", threads], stdout=sp.PIPE)
    write_to_log(mmseqs_search.stdout, logger)
    # creates the tsv output
    mmseqs_createtsv = sp.Popen(["mmseqs", "createtsv", os.path.join(phrog_db_dir, "phrogs_profile_db"), os.path.join(target_db_dir, "target_seqs"), os.path.join(mmseqs_dir, "results_mmseqs"), 
    os.path.join(out_dir,"mmseqs_results.tsv"), "--full-header", "--threads", threads], stdout=sp.PIPE)
    write_to_log(mmseqs_createtsv.stdout, logger)
    # remove the target dir when finished 
    sp.run(["rm", "-r", target_db_dir], check=True)


def convert_gff_to_gbk(filepath_in, out_dir, prefix):
    """
    Converts the gff to genbank
    :param filepath_in: input fasta file
    :param out_dir: output directory
    :param prefix: prefix
    :return:
    """
    gff_file = os.path.join(out_dir, prefix + ".gff")
    gbk_file = os.path.join(out_dir, prefix + ".gbk")
    with open(gbk_file, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(filepath_in, "fasta"))
        for record in GFF.parse(gff_file, fasta_handler):
            # instantiate record
            record.annotations["molecule_type"] = "DNA"
            record.annotations["date"] = datetime.today()
            record.annotations["topology"] = "linear"
            record.annotations["data_file_division"] = "VRL"
            # add features to the record
            for feature in record.features:
                # add translation only if CDS
                if feature.type == "CDS":
                    if feature.strand == 1:
                        feature.qualifiers.update({'translation': Seq.translate(record.seq[feature.location.start.position:feature.location.end.position], to_stop=True)})
                    else: # reverse strand -1 needs reverse compliment
                        feature.qualifiers.update({'translation': Seq.translate(record.seq[feature.location.start.position:feature.location.end.position].reverse_complement(), to_stop=True)})
            SeqIO.write(record, gbk_handler, "genbank")

def run_minced(filepath_in, out_dir, prefix, logger):
    """
    Runs MinCED 
    :param filepath_in: input fasta file
    :param out_dir: output directory
    :param logger: logger
    :params prefix: prefix
    :return:
    """
    print("Running MinCED.")
    logger.info("Running MinCED.")
    try:
        minced_fast = sp.Popen(["minced", filepath_in, os.path.join(out_dir, prefix + "_minced_spacers.txt") , os.path.join(out_dir, prefix + "_minced.gff")], stderr=sp.PIPE, stdout=sp.PIPE) 
        write_to_log(minced_fast.stderr, logger)
    except:
        sys.exit("Error with MinCED\n")  

def run_aragorn(filepath_in, out_dir, prefix, logger):
    """
    Runs run_aragorn 
    :param filepath_in: input fasta file
    :param out_dir: output directory
    :param logger: logger
    :params prefix: prefix
    :return:
    """
    print("Running Aragorn.")
    logger.info("Running Aragorn.")
    try:
        aragorn = sp.Popen(["aragorn", "-l", "-gcbact", "-w", "-o", os.path.join(out_dir, prefix + "_aragorn.txt"), "-m", filepath_in], stderr=sp.PIPE, stdout=sp.PIPE) 
        write_to_log(aragorn.stderr, logger)
    except:
        sys.exit("Error with Aragorn\n")  

def run_mmseqs_vfdb(db_dir, out_dir, threads, logger, gene_predictor):
    """
    Runs mmseqs2 on VFDB 
    :param db_dir: database path
    :param out_dir: output directory
    :param logger: logger
    :params threads: threads
    :param gene_predictor: phanotate or prodigal
    No evalue - settings as Enault et al. recommend 80% identify 40% coverage https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315482/ 
    :return:
    """
    print("Running mmseqs2 on vfdb.")
    logger.info("Running mmseqs2 on vfdb.")
    vfdb_db_dir = db_dir
    vfdb_dir = os.path.join(out_dir, "vfdb/")
    amino_acid_fasta = gene_predictor + "_aas_tmp.fasta"
    target_db_dir =  os.path.join(out_dir, "vfdb_target_dir/") 
    tmp_dir = os.path.join(out_dir, "vfdb_tmp_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input fasta
    vfdb_createdb = sp.Popen(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")], stdout=sp.PIPE)
    write_to_log(vfdb_createdb.stdout, logger)
    # runs the search
    vfdb_search = sp.Popen(["mmseqs", "search", "--min-seq-id", "0.8", "-c", "0.4", os.path.join(vfdb_db_dir, "vfdb"), os.path.join(target_db_dir, "target_seqs"), os.path.join(vfdb_dir, "results_vfdb"), tmp_dir, "-s", "8.5",
    "--threads", threads], stdout=sp.PIPE)
    write_to_log(vfdb_search.stdout, logger)
    # creates the tsv output
    vfdb_createtsv = sp.Popen(["mmseqs", "createtsv", os.path.join(vfdb_db_dir, "vfdb"), os.path.join(target_db_dir, "target_seqs"), os.path.join(vfdb_dir, "results_vfdb"), 
    os.path.join(out_dir,"vfdb_results.tsv"), "--full-header", "--threads", threads], stdout=sp.PIPE)
    write_to_log(vfdb_createtsv.stdout, logger)
    # remove the target dir when finished 
    sp.run(["rm", "-r", target_db_dir], check=True)


def run_mmseqs_card(db_dir, out_dir, threads, logger, gene_predictor):
    """
    Runs mmseqs2 on card 
    :param db_dir: database path
    :param out_dir: output directory
    :param logger: logger
    :params threads: threads
    :param gene_predictor: phanotate or prodigal
    No evalue - settings as Enault et al. recommend 80% identify 40% coverage https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5315482/ 
    :return:
    """
    print("Running mmseqs2 on CARD.")
    logger.info("Running mmseqs2 on CARD.")
    CARD_db_dir = db_dir
    CARD_dir = os.path.join(out_dir, "CARD/")
    amino_acid_fasta = gene_predictor + "_aas_tmp.fasta"
    target_db_dir =  os.path.join(out_dir, "CARD_target_dir/") 
    tmp_dir = os.path.join(out_dir, "CARD_tmp_dir/") 

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    CARD_createdb = sp.Popen(["mmseqs", "createdb", os.path.join(out_dir, amino_acid_fasta), os.path.join(target_db_dir, "target_seqs")], stdout=sp.PIPE)
    write_to_log(CARD_createdb.stdout, logger)
    # runs the seacrh
    CARD_search = sp.Popen(["mmseqs", "search", "--min-seq-id", "0.8", "-c", "0.4", os.path.join(CARD_db_dir, "CARD"), os.path.join(target_db_dir, "target_seqs"), os.path.join(CARD_dir, "results_CARD"), tmp_dir, "-s", "8.5",
    "--threads", threads], stdout=sp.PIPE)
    write_to_log(CARD_search.stdout, logger)
    # creates the tsv
    CARD_createtsv = sp.Popen(["mmseqs", "createtsv", os.path.join(CARD_db_dir, "CARD"), os.path.join(target_db_dir, "target_seqs"), os.path.join(CARD_dir, "results_CARD"), 
    os.path.join(out_dir,"CARD_results.tsv"), "--full-header", "--threads", threads], stdout=sp.PIPE)
    write_to_log(CARD_createtsv.stdout, logger)
    # remove the target dir when finished 
    sp.run(["rm", "-r", target_db_dir], check=True)



def reorient_terminase(filepath_in, out_dir, prefix, terminase_strand, terminase_start, logger):
    """
    re-orients phage to begin with large terminase subunit 
    :param filepath_in input genome fasta
    :param out_dir: output directory path
    :param prefix: prefix for pharokka
    :param terminase_strand: strandedness of the terminase large subunit. Is either 'pos' or 'neg'
    :param terminase_start: start coordinate of terminase large subunit. 
    :logger: logger
    """
    print("Input checked. \nReorienting input genome to begin with terminase large subunit.")
    logger.info("Input checked. \nReorienting input genome to begin with terminase large subunit.")


    # read in the fasta
    record = SeqIO.read(filepath_in, "fasta")

    # get length of the fasta
    length = len(record.seq)

    if int(terminase_start) > length or int(terminase_start) < 1:
        sys.exit("Error: terminase large subunit start coordinate specified is not within the provided genome length. Please check your input. \n")  

    # positive

    # reorient to start at the terminase  
    # pos
    if terminase_strand == "pos":
        start = record.seq[(int(terminase_start) - 1):length]
        end = record.seq[0:int(terminase_start)-1]
        total = start + end

    # neg
    if terminase_strand == "neg":
        record.seq = record.seq.reverse_complement()
        start = record.seq[(length - int(terminase_start)):length]
        end = record.seq[0:(length - int(terminase_start))]
        total = start + end

    # set sequence
    record.seq = total

    out_fasta = os.path.join(out_dir, prefix + '_genome_terminase_reoriented.fasta')

    SeqIO.write(record, out_fasta, "fasta")

def run_mash_sketch(filepath_in, out_dir, logger):
    """
    Runs mash sketch
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :return:
    """
    try:
        mash_sketch = sp.Popen(["mash", "sketch", filepath_in, "-o", os.path.join(out_dir, "input_mash_sketch.msh"), "-i"], stderr=sp.PIPE, stdout=sp.DEVNULL) 
        write_to_log(mash_sketch.stderr, logger)
    except:
        sys.exit("Error with mash sketch\n")  


def run_mash_dist( out_dir,db_dir, logger):
    """
    Runs mash
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :return:
    """
    mash_tsv = os.path.join(out_dir,"mash_out.tsv")
    outFile = open(mash_tsv, "w")
    try:
        mash_dist = sp.Popen(["mash", "dist", os.path.join(out_dir, "input_mash_sketch.msh"), os.path.join(db_dir, "5Jan2023_genomes.fa.msh"), "-d", "0.2", "-i" ], stdout=outFile, stderr=sp.PIPE) 
        write_to_log(mash_dist.stderr, logger)
    except:
        sys.exit("Error with mash dist\n")  


