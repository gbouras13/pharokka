import multiprocessing.pool
import os
import subprocess as sp
from datetime import datetime

import pandas as pd
import pyrodigal
import pyrodigal_gv
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from external_tools import ExternalTool
from loguru import logger
from util import remove_directory


def run_pyrodigal_gv(filepath_in, out_dir, threads):
    """
    Gets CDS using pyrodigal_gv
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :param meta Boolean - metagenomic mode flag
    :param coding_table coding table for prodigal (default 11)
    :return:
    """

    # true
    orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)

    def _find_genes(record):
        genes = orf_finder.find_genes(str(record.seq))
        return (record.id, genes)

    with multiprocessing.pool.ThreadPool(threads) as pool:
        with open(os.path.join(out_dir, "prodigal-gv_out.gff"), "w") as gff:
            with open(os.path.join(out_dir, "prodigal-gv_out_tmp.fasta"), "w") as dst:
                with open(
                    os.path.join(out_dir, "prodigal-gv_out_aas_tmp.fasta"), "w"
                ) as aa_fasta:
                    records = SeqIO.parse(filepath_in, "fasta")
                    for record_id, genes in pool.imap(_find_genes, records):
                        genes.write_gff(
                            gff, sequence_id=record_id, include_translation_table=True
                        )
                        genes.write_genes(dst, sequence_id=record_id)
                        # need to write the translation
                        genes.write_translations(aa_fasta, sequence_id=record_id)


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

    for i in range(1, num_fastas + 1):
        in_file = "input_subprocess" + str(i) + ".fasta"
        out_file = "phanotate_out_tmp" + str(i) + ".fasta"
        filepath_in = os.path.join(phanotate_tmp_dir, in_file)
        cmd = (
            "phanotate.py "
            + filepath_in
            + " -o "
            + os.path.join(phanotate_tmp_dir, out_file)
            + " -f fasta"
        )
        commands.append(cmd)

    n = int(threads)  # the number of parallel processes you want

    for j in range(max(int(len(commands) / n) + 1, 1)):
        procs = [
            sp.Popen(i, shell=True)
            for i in commands[j * n : min((j + 1) * n, len(commands))]
        ]
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

    for i in range(1, num_fastas + 1):
        in_file = "input_subprocess" + str(i) + ".fasta"
        out_file = "phanotate_out_tmp" + str(i) + ".txt"
        filepath_in = os.path.join(phanotate_tmp_dir, in_file)
        cmd = (
            "phanotate.py "
            + filepath_in
            + " -o "
            + os.path.join(phanotate_tmp_dir, out_file)
            + " -f tabular"
        )
        commands.append(cmd)

    n = int(threads)  # the number of parallel processes you want
    for j in range(max(int(len(commands) / n) + 1, 1)):
        procs = [
            sp.Popen(i, shell=True)
            for i in commands[j * n : min((j + 1) * n, len(commands))]
        ]
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
    for i in range(1, int(num_fastas) + 1):
        out_tsv = "phanotate_out_tmp" + str(i) + ".txt"
        tsvs.append(os.path.join(phanotate_tmp_dir, out_tsv))

    with open(os.path.join(out_dir, "phanotate_out.txt"), "w") as outfile:
        for fname in tsvs:
            with open(fname) as infile:
                outfile.write(infile.read())

    fastas = []
    for i in range(1, int(num_fastas) + 1):
        out_fasta = "phanotate_out_tmp" + str(i) + ".fasta"
        fastas.append(os.path.join(phanotate_tmp_dir, out_fasta))

    with open(os.path.join(out_dir, "phanotate_out_tmp.fasta"), "w") as outfile:
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

    for i in range(1, num_fastas + 1):
        in_file = "input_subprocess" + str(i) + ".fasta"
        out_file = "trnascan_tmp" + str(i) + ".gff"
        filepath_in = os.path.join(input_tmp_dir, in_file)
        filepath_out = os.path.join(input_tmp_dir, out_file)
        cmd = "tRNAscan-SE " + filepath_in + " --thread 1 -G -Q -j " + filepath_out
        commands.append(cmd)

    n = int(threads)  # the number of parallel processes you want

    for j in range(max(int(len(commands) / n) + 1, 1)):
        procs = [
            sp.Popen(i, shell=True, stderr=sp.PIPE, stdout=sp.DEVNULL)
            for i in commands[j * n : min((j + 1) * n, len(commands))]
        ]
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
    for i in range(1, int(num_fastas) + 1):
        out_gff = "trnascan_tmp" + str(i) + ".gff"
        gffs.append(os.path.join(input_tmp_dir, out_gff))

    with open(os.path.join(out_dir, "trnascan_out.gff"), "w") as outfile:
        for fname in gffs:
            with open(fname) as infile:
                outfile.write(infile.read())


##### single contig mode ######


def run_phanotate(filepath_in, out_dir, logdir):
    """
    Runs phanotate
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logdir logdir
    :return:
    """

    out_fasta = os.path.join(out_dir, "phanotate_out_tmp.fasta")
    out_tab = os.path.join(out_dir, "phanotate_out.txt")

    phan_fast = ExternalTool(
        tool="phanotate.py",
        input=f"{filepath_in}",
        output=f"-o {out_fasta}",
        params=f"-f fasta",
        logdir=logdir,
        outfile="",
    )

    phan_txt = ExternalTool(
        tool="phanotate.py",
        input=f"{filepath_in}",
        output=f"-o {out_tab}",
        params=f"-f tabular",
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(phan_fast)
        ExternalTool.run_tool(phan_txt)

    except:
        logger.error("Error with Phanotate\n")


def run_pyrodigal(filepath_in, out_dir, meta, coding_table, threads):
    """
    Gets CDS using pyrodigal
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :param meta Boolean - metagenomic mode flag
    :param coding_table coding table for prodigal (default 11)
    :param threads: threads
    :return:
    """

    prodigal_metamode = False
    if meta == True:
        prodigal_metamode = True
        logger.info("Prodigal Meta Mode Enabled")

    #######################################################
    # if under 20000, pyrodigal will only work in meta mode
    # https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#plasmids-phages-viruses-and-other-short-sequences
    # https://github.com/hyattpd/Prodigal/issues/51
    # so make sure of this
    #######################################################

    # get total length of input
    total_length = 0

    with open(filepath_in, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total_length += len(record.seq)

    # if the length is 100000 or under, use meta mode by default
    if total_length < 100001:
        orf_finder = pyrodigal.GeneFinder(meta=True)
    # otherwise train it
    # recommend pyrodigal-gv anyway
    else:
        # for training if you want different coding table
        seqs = [bytes(record.seq) for record in SeqIO.parse(filepath_in, "fasta")]
        record = SeqIO.parse(filepath_in, "fasta")
        orf_finder = pyrodigal.GeneFinder(meta=prodigal_metamode)

        # make coding table possible if false
        if prodigal_metamode == False:
            orf_finder.train(*seqs, translation_table=int(coding_table))

    # define for the multithreadpool
    def _find_genes(record):
        genes = orf_finder.find_genes(str(record.seq))
        return (record.id, genes)

    with multiprocessing.pool.ThreadPool(threads) as pool:
        with open(os.path.join(out_dir, "prodigal_out.gff"), "w") as gff:
            with open(os.path.join(out_dir, "prodigal_out_tmp.fasta"), "w") as dst:
                with open(
                    os.path.join(out_dir, "prodigal_out_aas_tmp.fasta"), "w"
                ) as aa_fasta:
                    records = SeqIO.parse(filepath_in, "fasta")
                    for record_id, genes in pool.imap(_find_genes, records):
                        genes.write_gff(gff, sequence_id=record_id)
                        genes.write_genes(dst, sequence_id=record_id)
                        # need to write the translation
                        genes.write_translations(aa_fasta, sequence_id=record_id)


def tidy_phanotate_output(out_dir):
    """
    Tidies phanotate output
    :param out_dir: output directory
    :return: phan_df pandas dataframe
    """
    phan_file = os.path.join(out_dir, "phanotate_out.txt")
    col_list = ["start", "stop", "frame", "contig", "score"]
    dtype_dict = {
        "start": int,
        "stop": int,
        "frame": str,
        "contig": str,
        "score": float,
    }

    phan_df = pd.read_csv(
        phan_file,
        delimiter="\t",
        index_col=False,
        names=col_list,
        skiprows=2, # to skip the headers
        dtype=dtype_dict,
        comment="#"
    )
    # get rid of the headers and reset the index
    phan_df = phan_df[phan_df["start"] != "#id:"]
    phan_df = phan_df[phan_df["start"] != "#START"].reset_index(drop=True)
    phan_df["gene"] = (
        phan_df["contig"].astype(str)
        + phan_df.index.astype(str)
        + " "
        + phan_df["start"].astype(str)
        + "_"
        + phan_df["stop"].astype(str)
    )
    phan_df.to_csv(
        os.path.join(out_dir, "cleaned_phanotate.tsv"), sep="\t", index=False
    )
    return phan_df


def tidy_prodigal_output(out_dir, gv_flag):
    """
    Tidies prodigal output
    :param out_dir: output directory
    :param gv_flag: if prodigal-gv, then True
    :return: prod_filt_df pandas dataframe
    """
    if gv_flag is True:
        prefix = "prodigal-gv"
    else:
        prefix = "prodigal"

    prod_file = os.path.join(out_dir, f"{prefix}_out.gff")
    col_list = [
        "contig",
        "prod",
        "orf",
        "start",
        "stop",
        "score",
        "frame",
        "phase",
        "description",
    ]
    dtype_dict = {
        "contig": str,
        "prod": str,
        "orf": str,
        "start": int,
        "stop": int,
        "score": float,
        "frame": str,
        "phase": str,
        "description": str,
    }

    prod_df = pd.read_csv(
        prod_file,
        delimiter="\t",
        index_col=False,
        names=col_list,
        dtype=dtype_dict,
        comment="#",  # to skip the headers
    )

    # meta mode brings in some Nas so remove them
    # need to reset index!!!! and drop, or else will cause rubbish results for metagenomics
    prod_df = prod_df.dropna().reset_index(drop=True)

    prod_filt_df = prod_df[["start", "stop", "frame", "contig", "score"]]

    # convert start stop to int
    prod_filt_df["start"] = prod_filt_df["start"].astype("int")
    prod_filt_df["stop"] = prod_filt_df["stop"].astype("int")
    # rearrange start and stop so that for negative strand, the stop is before start (like phanotate_out)
    cols = ["start", "stop"]
    # indices where start is greater than stop
    ixs = prod_filt_df["frame"] == "-"
    # Where ixs is True, values are swapped
    prod_filt_df.loc[ixs, cols] = (
        prod_filt_df.loc[ixs, cols].reindex(columns=cols[::-1]).values
    )
    prod_filt_df["gene"] = (
        prod_filt_df["contig"].astype(str)
        + prod_filt_df.index.astype(str)
        + " "
        + prod_filt_df["start"].astype(str)
        + "_"
        + prod_filt_df["stop"].astype(str)
    )
    prod_filt_df.to_csv(
        os.path.join(out_dir, f"cleaned_{prefix}.tsv"), sep="\t", index=False
    )
    return prod_filt_df


def tidy_genbank_output(out_dir, genbank_file, coding_table):
    """
    Tidies phanotate output
    :param out_dir: output directory
    :return: phan_df pandas dataframe
    """

    # Create lists to store CDS information
    starts = []
    stops = []
    frames = []
    contigs = []

    # output nuc
    fasta_nucl_tmp = "genbank_out_tmp.fasta"
    # output aa file
    fasta_output_aas_tmp = "genbank_aas_tmp.fasta"

    # make dataframe
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                frame = feature.location.strand
                if frame == 1:  # pos
                    frame = "+"
                    start = feature.location.start
                    stop = feature.location.end
                else:  # neg
                    frame = "-"
                    start = feature.location.end
                    stop = feature.location.start

                contig = record.id
                starts.append(start)
                stops.append(stop)
                frames.append(frame)
                contigs.append(contig)

    # Create a pandas DataFrame
    data = {"start": starts, "stop": stops, "frame": frames, "contig": contigs}

    gen_df = pd.DataFrame(data)
    # add fake score
    gen_df["score"] = "No_score"

    # get the gene
    gen_df["gene"] = (
        gen_df["contig"].astype(str)
        + gen_df.index.astype(str)
        + " "
        + gen_df["start"].astype(str)
        + "_"
        + gen_df["stop"].astype(str)
    )

    gen_df.to_csv(os.path.join(out_dir, "cleaned_genbank.tsv"), sep="\t", index=False)

    # get temporary AA output
    with open(os.path.join(out_dir, fasta_nucl_tmp), "w") as nt_fa:
        with open(os.path.join(out_dir, fasta_output_aas_tmp), "w") as aa_fa:
            i = 0
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        header = str(gen_df["contig"].iloc[i]) + str(i)
                        description = (
                            str(gen_df["start"].iloc[i])
                            + "_"
                            + str(gen_df["stop"].iloc[i])
                        )

                        # Extract NT and CDS sequence
                        cds_nt_seq = feature.location.extract(record.seq)
                        cds_amino_acid_seq = cds_nt_seq.translate(
                            to_stop=True, table=coding_table
                        )

                        # save to file
                        nt_record = aa_record = SeqRecord(
                            cds_nt_seq,
                            id=header,
                            description=description,
                        )

                        aa_record = aa_record = SeqRecord(
                            cds_amino_acid_seq,
                            id=header,
                            description=description,
                        )

                        SeqIO.write(aa_record, aa_fa, "fasta")
                        SeqIO.write(nt_record, nt_fa, "fasta")
                        i += 1

    return gen_df


def translate_fastas(out_dir, gene_predictor, coding_table, genbank_file):
    """
    Translates input CDSs to amino acids. For genbank, it will just output the CDS coordinates in the file
    :param out_dir: output directory
    :param gene_predictor: phanotate or prodigal
    :param genbank_file: genbank file if gene_predictor is genbank
    :return:
    """
    if gene_predictor == "phanotate":
        clean_df = tidy_phanotate_output(out_dir)
    elif gene_predictor == "prodigal":
        clean_df = tidy_prodigal_output(out_dir, False)  # gv_flag is false
    elif gene_predictor == "prodigal-gv":
        clean_df = tidy_prodigal_output(out_dir, True)  # gv_flag is true
    elif gene_predictor == "genbank":
        clean_df = tidy_genbank_output(out_dir, genbank_file, coding_table)

    fasta_output_aas_tmp = gene_predictor + "_aas_tmp.fasta"

    if gene_predictor == "phanotate":
        # read the nucl fasta
        fasta_input_tmp = gene_predictor + "_out_tmp.fasta"
        # translate for temporary AA output
        with open(os.path.join(out_dir, fasta_output_aas_tmp), "w") as aa_fa:
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(out_dir, fasta_input_tmp), "fasta"
            ):
                dna_header = str(clean_df["contig"].iloc[i]) + str(i)
                dna_description = (
                    str(clean_df["start"].iloc[i]) + "_" + str(clean_df["stop"].iloc[i])
                )
                aa_record = SeqRecord(
                    dna_record.seq.translate(to_stop=True, table=coding_table),
                    id=dna_header,
                    description=dna_description,
                )
                SeqIO.write(aa_record, aa_fa, "fasta")
                i += 1
    elif gene_predictor == "prodigal-gv" or gene_predictor == "prodigal":
        # read in the AA file instead and parse that to clean the header
        fasta_input_tmp = gene_predictor + "_out_aas_tmp.fasta"
        with open(os.path.join(out_dir, fasta_output_aas_tmp), "w") as aa_fa:
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(out_dir, fasta_input_tmp), "fasta"
            ):
                dna_header = str(clean_df["contig"].iloc[i]) + str(i)
                dna_description = (
                    str(clean_df["start"].iloc[i]) + "_" + str(clean_df["stop"].iloc[i])
                )
                aa_record = SeqRecord(
                    dna_record.seq,
                    id=dna_header,
                    description=dna_description,
                )
                SeqIO.write(aa_record, aa_fa, "fasta")
                i += 1
    # for genbank do nothing


def run_trna_scan(filepath_in, threads, out_dir, logdir):
    """
    Runs trna scan
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :return:
    """

    out_gff = os.path.join(out_dir, "trnascan_out.gff")

    trna = ExternalTool(
        tool="tRNAscan-SE",
        input=f"{filepath_in}",
        output=f"{out_gff}",
        params=f"--thread {threads} -G -Q -j",
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(trna)
    except:
        logger.error("Error with tRNAscan-SE")
        return 0


def run_mmseqs(db_dir, out_dir, threads, logdir, gene_predictor, evalue, db_name):
    """
    Runs mmseqs2 on phrogs
    :param db_dir: database path
    :param out_dir: output directory
    :param logger: logger
    :params threads: threads
    :param gene_predictor: phanotate or prodigal
    :param evalue: evalue for mmseqs2
    :param db_name: str one of 'PHROG', 'VFDB' or 'CARD'
    :return:
    """

    logger.info(f"Running MMseqs2 on {db_name} Database.")

    # declare files
    amino_acid_fasta = gene_predictor + "_aas_tmp.fasta"

    # define the outputs
    if db_name == "PHROG":
        mmseqs_dir = os.path.join(out_dir, "mmseqs/")
        target_db_dir = os.path.join(out_dir, "target_dir/")
        tmp_dir = os.path.join(out_dir, "tmp_dir/")
        profile_db = os.path.join(db_dir, "phrogs_profile_db")
        mmseqs_result_tsv = os.path.join(out_dir, "mmseqs_results.tsv")
    elif db_name == "VFDB":
        mmseqs_dir = os.path.join(out_dir, "VFDB/")
        target_db_dir = os.path.join(out_dir, "VFDB_target_dir/")
        tmp_dir = os.path.join(out_dir, "VFDB_dir/")
        profile_db = os.path.join(db_dir, "vfdb")
        mmseqs_result_tsv = os.path.join(out_dir, "vfdb_results.tsv")
    elif db_name == "CARD":
        mmseqs_dir = os.path.join(out_dir, "CARD/")
        target_db_dir = os.path.join(out_dir, "CARD_target_dir/")
        tmp_dir = os.path.join(out_dir, "CARD_dir/")
        profile_db = os.path.join(db_dir, "CARD")
        mmseqs_result_tsv = os.path.join(out_dir, "CARD_results.tsv")

    input_aa_fasta = os.path.join(out_dir, amino_acid_fasta)
    target_seqs = os.path.join(target_db_dir, "target_seqs")

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    mmseqs_createdb = ExternalTool(
        tool="mmseqs createdb",
        input=f"",
        output=f"{target_seqs}",
        params=f"{input_aa_fasta}",  # param goes before output and mmseqs2 required order
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(mmseqs_createdb)

    # runs the mmseqs seacrh

    result_mmseqs = os.path.join(mmseqs_dir, "results_mmseqs")

    if db_name == "PHROG":
        mmseqs_search = ExternalTool(
            tool="mmseqs search",
            input=f"",
            output=f"{tmp_dir} -s 8.5 --threads {threads}",
            params=f"-e {evalue} {profile_db} {target_seqs} {result_mmseqs}",  # param goes before output and mmseqs2 required order
            logdir=logdir,
            outfile="",
        )
    else:  # if it is vfdb or card search with cutoffs instead of evalue
        mmseqs_search = ExternalTool(
            tool="mmseqs search",
            input=f"",
            output=f"{tmp_dir} -s 8.5 --threads {threads}",
            params=f"--min-seq-id 0.8 -c 0.4 {profile_db} {target_seqs} {result_mmseqs}",  # param goes before output and mmseqs2 required order
            logdir=logdir,
            outfile="",
        )

    ExternalTool.run_tool(mmseqs_search)

    # creates the output tsv
    mmseqs_createtsv = ExternalTool(
        tool="mmseqs createtsv",
        input=f"",
        output=f"{mmseqs_result_tsv} --full-header --threads {threads}",
        params=f"{profile_db} {target_seqs} {result_mmseqs} ",  # param goes before output and mmseqs2 required order
        logdir=logdir,
        outfile="",
    )

    ExternalTool.run_tool(mmseqs_createtsv)

    # remove the target dir when finished
    remove_directory(target_db_dir)


def convert_gff_to_gbk(filepath_in, input_dir, out_dir, prefix, prot_seq_df):
    """
    Converts the gff to genbank
    :param filepath_in: input fasta file
    :param input_dir: input directory of the gff. same as output_dir for the overall gff in normal mode, differeny for meta mode
    :param out_dir: output directory of the gbk
    :param prefix: prefix
    :param prefix: prot_seq_df from pharok object with gene name + protein sequence for all genes (from create_gff()).
    :return:
    """

    gff_file = os.path.join(input_dir, f"{prefix}.gff")
    gbk_file = os.path.join(out_dir, f"{prefix}.gbk")

    prot_seq_df["contig"] = prot_seq_df["contig"].astype(str)

    with open(gbk_file, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(filepath_in, "fasta"))
        for record in GFF.parse(gff_file, fasta_handler):
            # sequence in each contig (record)
            record.id = str(record.id)
            subset_seqs_df = prot_seq_df.loc[prot_seq_df["contig"] == record.id]
            # get all the seqs in the contigs - and drop the index to reset for 0 indexed loop
            subset_seqs = subset_seqs_df["sequence"].reset_index(drop=True)
            # start the loop
            i = 0

            # instantiate record
            record.annotations["molecule_type"] = "DNA"
            record.annotations["date"] = datetime.today()
            record.annotations["topology"] = "linear"
            record.annotations[
                "data_file_division"
            ] = "PHG"  # https://github.com/RyanCook94/inphared/issues/22
            # add features to the record
            for feature in record.features:
                # add translation only if CDS
                if feature.type == "CDS":
                    # aa = prot_records[i].seq
                    if feature.strand == 1:
                        feature.qualifiers.update(
                            {"translation": subset_seqs[i]}  # from the aa seq
                        )
                    else:  # reverse strand -1 needs reverse compliment
                        feature.qualifiers.update(
                            {"translation": subset_seqs[i]}  # from the aa seq
                        )
                    i += 1
            SeqIO.write(record, gbk_handler, "genbank")


def run_minced(filepath_in, out_dir, prefix, minced_args, logdir):
    """
    Runs MinCED
    :param filepath_in: input fasta file
    :param out_dir: output directory
    :param logger: logger
    :param minced_args: str with extra arguments to pass to MINced
    :params prefix: prefix
    :return:
    """

    logger.info("Running MinCED.")

    output_spacers = os.path.join(out_dir, prefix + "_minced_spacers.txt")
    output_gff = os.path.join(out_dir, prefix + "_minced.gff")

    if minced_args != "":
        minced_args = f"-{minced_args}"

    minced_fast = ExternalTool(
        tool="minced",
        input=f"",
        output=f" {output_spacers} {output_gff}",
        params=f" {minced_args} {filepath_in}",  # need strange order for minced params go first
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(minced_fast)
    except:
        logger.error("Error with MinCED\n")


def run_aragorn(filepath_in, out_dir, prefix, logdir):
    """
    Runs run_aragorn
    :param filepath_in: input fasta file
    :param out_dir: output directory
    :param logdir: logdir
    :params prefix: prefix
    :return:
    """
    logger.info("Running Aragorn.")

    aragorn_out_file = os.path.join(out_dir, prefix + "_aragorn.txt")
    aragorn = ExternalTool(
        tool="aragorn",
        input=f"{filepath_in}",
        output=f"-o {aragorn_out_file}",
        params=f"-l -gcbact -w -m",
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(aragorn)
    except:
        logger.error("Error with Aragorn\n")


def reorient_terminase(filepath_in, out_dir, prefix, terminase_strand, terminase_start):
    """
    re-orients phage to begin with large terminase subunit
    :param filepath_in input genome fasta
    :param out_dir: output directory path
    :param prefix: prefix for pharokka
    :param terminase_strand: strandedness of the terminase large subunit. Is either 'pos' or 'neg'
    :param terminase_start: start coordinate of terminase large subunit.
    :logger: logger
    """

    logger.info(
        "Input checked. \nReorienting input genome to begin with terminase large subunit."
    )

    # read in the fasta
    record = SeqIO.read(filepath_in, "fasta")

    # get length of the fasta
    length = len(record.seq)

    if int(terminase_start) > length or int(terminase_start) < 1:
        logger.error(
            "Error: terminase large subunit start coordinate specified is not within the provided genome length. Please check your input. \n"
        )

    # positive

    # reorient to start at the terminase
    # pos
    if terminase_strand == "pos":
        start = record.seq[(int(terminase_start) - 1) : length]
        end = record.seq[0 : int(terminase_start) - 1]
        total = start + end

    # neg
    if terminase_strand == "neg":
        record.seq = record.seq.reverse_complement()
        start = record.seq[(length - int(terminase_start)) : length]
        end = record.seq[0 : (length - int(terminase_start))]
        total = start + end

    # set sequence
    record.seq = total

    out_fasta = os.path.join(out_dir, prefix + "_genome_terminase_reoriented.fasta")

    SeqIO.write(record, out_fasta, "fasta")


def run_mash_sketch(filepath_in, out_dir, logdir):
    """
    Runs mash sketch
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param logger logger
    :return:
    """

    mash_sketch_out_file = os.path.join(out_dir, "input_mash_sketch.msh")

    mash_sketch = ExternalTool(
        tool="mash sketch",
        input=f"",
        output=f"-o {mash_sketch_out_file} -i",
        params=f"{filepath_in}",
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(mash_sketch)

    except:
        logger.error("Error with mash sketch\n")


def run_mash_dist(out_dir, db_dir, mash_distance, logdir):
    """
    Runs mash
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param mash_distance: mash distance - float
    :param logger logger
    :return:
    """

    mash_sketch = os.path.join(out_dir, "input_mash_sketch.msh")
    phrog_sketch = os.path.join(db_dir, "1Aug2023_genomes.fa.msh")
    mash_tsv = os.path.join(out_dir, "mash_out.tsv")

    mash_dist = ExternalTool(
        tool="mash",
        input="",
        output="",
        params=f" dist  {mash_sketch} {phrog_sketch} -d {mash_distance} -i ",
        logdir=logdir,
        outfile=mash_tsv,
    )

    # need to write to stdout

    try:
        ExternalTool.run_tool(mash_dist, to_stdout=True)

    except:
        logger.error("Error with mash dist\n")


def run_dnaapler(filepath_in, contig_count, out_dir, threads, logdir):
    """
    Runs dnaapler
    :param filepath_in: input fasta file
    :param contig_count: int numbher of contigs
    :param out_dir: output directory
    :param threads: number of threads
    :params prefix: prefix
    :return: dnaapler_success - whether dnaapler worked
    """

    logger.info(
        "Running Dnaapler to rerorient all contigs to begin with the terminase large subunit."
    )

    dnaapler_output_dir = os.path.join(out_dir, "dnaapler")

    if contig_count == 1:
        dnaapler = ExternalTool(
            tool="dnaapler phage",
            input=f"-i {filepath_in}",
            output=f"-o {dnaapler_output_dir} -t {threads}",
            params=f"",  # need strange order for minced params go first
            logdir=logdir,
            outfile="",
        )
    else:
        dnaapler = ExternalTool(
            tool="dnaapler all",
            input=f"-i {filepath_in}",
            output=f"-o {dnaapler_output_dir} -t {threads}",
            params=f"",  # need strange order for minced params go first
            logdir=logdir,
            outfile="",
        )

    dnaapler_success = True
    try:
        ExternalTool.run_tool(dnaapler)
    except:
        logger.warning("Dnaapler failed to find the large terminase subunit.")
        logger.warning(
            "For reorientation, please run pharokka again without --dnaapler and use --terminase instead."
        )
        logger.warning("Pharokka will now continue without reorientation.")
        dnaapler_success = False

    return dnaapler_success


# #### skani

# # no need to sketch
# # takes like 10 secs to run this

# # pharokka_v1.4.0_databases % skani dist -q ../pharokka/tests/test_data/overall/Meta_example/combined_meta.fasta -r 1Jan2024_genomes.fa  --qi --ri --slow -m 200 -n 1  > out


# def run_skani(input_fasta, out_dir, db_dir, logdir):
#     """
#     Runs mash
#     :param filepath_in: input filepath
#     :param out_dir: output directory
#     :param mash_distance: mash distance - float
#     :param logger logger
#     :return:
#     """

#     inphared_fasta = os.path.join(db_dir, "1Jan2024_genomes.fa.gz")
#     skani_tsv = os.path.join(out_dir, "skani_out.tsv")

#     mash_dist = ExternalTool(
#         tool="skani",
#         input="",
#         output="",
#         params=f" dist  -q {input_fasta} -r {inphared_fasta} --qi --ri --slow -m 200 -n 1 ",
#         logdir=logdir,
#         outfile=skani_tsv,
#     )

#     # need to write to stdout

#     try:
#         ExternalTool.run_tool(mash_dist, to_stdout=True)

#     except:
#         logger.error("Error with mash dist\n")
