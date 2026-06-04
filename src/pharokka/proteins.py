"""
Holds functions for pharokka proteins command
"""

import argparse
import collections
import os
import sys
import time
from argparse import RawTextHelpFormatter
from pathlib import Path

import polars as pl
import pyhmmer
from Bio import SeqIO
from loguru import logger

from .databases import check_db_installation
from .external_tools import ExternalTool
from .input_commands import check_dependencies, instantiate_dirs, validate_fasta, validate_threads
from .post_processing import process_card_results, process_pyhmmer_results, process_vfdb_results
from .util import get_contig_headers, get_version, remove_directory, remove_file

Result = collections.namedtuple("Result", ["protein", "phrog", "bitscore", "evalue"])


def get_input_proteins():
    parser = argparse.ArgumentParser(
        description="pharokka proteins: fast phage protein annotation script",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        action="store",
        help="Input proteins file in amino acid fasta format.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Directory to write the output to.",
        default=os.path.join(os.getcwd(), "output/"),
    )
    parser.add_argument(
        "-d",
        "--database",
        action="store",
        help="Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.",
        default="Default",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads. Defaults to 1.",
        action="store",
        default=str(1),
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        help="Prefix for output files. This is not required.",
        default="Default",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        help="E-value threshold for MMseqs2 database PHROGs, VFDB and CARD and pyhmmer PHROGs database search. Defaults to 1E-05.",
        action="store",
        default="1E-05",
    )
    parser.add_argument(
        "--hmm_only",
        help="Runs pyhmmer (HMMs) with PHROGs only, not MMseqs2 with PHROGs, CARD or VFDB.",
        action="store_true",
    )
    parser.add_argument(
        "--mmseqs2_only",
        help="Runs MMseqs2 with PHROGs, CARD and VFDB only (same as Pharokka v1.3.2 and prior).",
        action="store_true",
    )
    parser.add_argument(
        "--reverse_mmseqs2",
        help="MMseqs2 database as target not query.",
        action="store_true",
    )
    parser.add_argument(
        "--sensitivity",
        help="MMseqs2 sensitivity.",
        default=8.5,
        type=float,
    )
    parser.add_argument(
        "-V",
        "--version",
        help="Print pharokka Version",
        action="version",
        version=get_version(),
    )
    parser.add_argument(
        "--citation", help="Print pharokka Citation", action="store_true"
    )
    args = parser.parse_args()

    return args


def run_mmseqs_proteins(input_fasta, db_dir, out_dir, threads, logdir, evalue, reverse_mmseqs2, sensitivity, db_name):
    """
    Runs mmseqs2 for pharokka proteins
    """

    logger.info(f"Running MMseqs2 on {db_name} Database.")

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
        tmp_dir = os.path.join(out_dir, "VFDB_dir/")
        profile_db = os.path.join(db_dir, "CARD")
        mmseqs_result_tsv = os.path.join(out_dir, "CARD_results.tsv")

    target_seqs = os.path.join(target_db_dir, "target_seqs")

    if not os.path.isdir(target_db_dir):
        os.mkdir(target_db_dir)

    mmseqs_createdb = ExternalTool(
        tool="mmseqs createdb",
        input=f"",
        output=f"{target_seqs}",
        params=f"{input_fasta}",
        logdir=logdir,
        outfile="",
    )
    ExternalTool.run_tool(mmseqs_createdb)

    result_mmseqs = os.path.join(mmseqs_dir, "results_mmseqs")

    if db_name == "PHROG":
        if reverse_mmseqs2:
            params_list = f"-e {evalue} {target_seqs} {profile_db} {result_mmseqs}"
        else:
            params_list = f"-e {evalue} {profile_db} {target_seqs} {result_mmseqs}"

        mmseqs_search = ExternalTool(
            tool="mmseqs search",
            input=f"",
            output=f"{tmp_dir} -s {sensitivity} --threads {threads}",
            params=params_list,
            logdir=logdir,
            outfile="",
        )
    else:
        if reverse_mmseqs2:
            params_list = f"--min-seq-id 0.8 -c 0.4 {target_seqs} {profile_db} {result_mmseqs}"
        else:
            params_list = f"--min-seq-id 0.8 -c 0.4 {profile_db} {target_seqs} {result_mmseqs}"

        mmseqs_search = ExternalTool(
            tool="mmseqs search",
            input=f"",
            output=f"{tmp_dir} -s {sensitivity} --threads {threads}",
            params=params_list,
            logdir=logdir,
            outfile="",
        )

    ExternalTool.run_tool(mmseqs_search)

    if reverse_mmseqs2:
        params_list = f"{target_seqs} {profile_db} {result_mmseqs}"
    else:
        params_list = f"{profile_db} {target_seqs} {result_mmseqs}"

    mmseqs_createtsv = ExternalTool(
        tool="mmseqs createtsv",
        input=f"",
        output=f"{mmseqs_result_tsv} --full-header --threads {threads}",
        params=params_list,
        logdir=logdir,
        outfile="",
    )
    ExternalTool.run_tool(mmseqs_createtsv)

    remove_directory(target_db_dir)


def run_pyhmmer_proteins(input_fasta, db_dir, threads, evalue):
    """
    Runs pyhmmer on phrogs
    """

    Result = collections.namedtuple(
        "Result", ["protein", "phrog", "bitscore", "evalue"]
    )

    results = []
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(os.path.join(db_dir, "all_phrogs.h3m"), alphabet=alphabet) as hmms:
        with pyhmmer.easel.SequenceFile(
            input_fasta, digital=True, alphabet=alphabet
        ) as seqs:
            for hits in pyhmmer.hmmer.hmmscan(
                seqs, hmms, cpus=int(threads), E=float(evalue)
            ):
                protein = hits.query.name
                for hit in hits:
                    if hit.included:
                        results.append(
                            Result(protein, hit.name, hit.score, hit.evalue)
                        )

    best_results = {}
    for result in results:
        if result.protein not in best_results or result.bitscore > best_results[result.protein].bitscore:
            best_results[result.protein] = result

    return best_results


class Pharok_Prot:
    """pharokka proteins Output Class"""

    def __init__(
        self,
        out_dir: Path = "output/",
        db_dir: Path = "db_dir/",
        prefix: str = "pharokka",
        input_fasta: Path = "input.fasta",
        pyhmmer_results_dict: dict = None,
        tophits_df=None,
        vfdb_results=None,
        card_results=None,
        length_df=None,
        mmseqs_flag: bool = True,
        hmm_flag: bool = True,
        reverse_mmseqs2: bool = False,
    ) -> None:
        """
        Parameters
        --------
        out_dir: str,
            output directory
        db_dir: str,
            database directory
        prefix: str, prefix for output
            prefix
        input_fasta: Path
            input FASTA file with phage contigs
        pyhmmer_results_dict: dict
            dictionary of Result tuples from pyhmmer.
        tophits_df: pl.DataFrame, required
            merged dataframe output
        vfdb_results: pl.DataFrame,
            vfdb dataframe output
        card_results: pl.DataFrame,
            CARD dataframe output
        length_df: pl.DataFrame,
            dataframe with lengths for each input contig
        mmseqs_flag: bool
            whether MMseqs2 was run
        hmm_flag: bool
            whether HMM was run
        reverse_mmseqs2: bool
            whether to run MMseqs using database as target
        """
        self.out_dir = out_dir
        self.db_dir = db_dir
        self.prefix = prefix
        self.input_fasta = input_fasta
        self.pyhmmer_results_dict = pyhmmer_results_dict if pyhmmer_results_dict is not None else {}
        self.tophits_df = tophits_df if tophits_df is not None else pl.DataFrame()
        self.vfdb_results = vfdb_results if vfdb_results is not None else pl.DataFrame()
        self.card_results = card_results if card_results is not None else pl.DataFrame()
        self.length_df = length_df if length_df is not None else pl.DataFrame()
        self.mmseqs_flag = mmseqs_flag
        self.hmm_flag = hmm_flag
        self.reverse_mmseqs2 = reverse_mmseqs2

    def process_dataframes(self):
        """
        Processes and combines PHROGS, CARD and VFDB mmseqs output
        """

        ids = get_contig_headers(self.input_fasta)  # returns polars Series (seq descriptions)
        all_genes_df = pl.DataFrame({"gene": ids.to_list()})

        if self.mmseqs_flag is True:
            mmseqs_file = os.path.join(self.out_dir, "mmseqs_results.tsv")
            logger.info("Processing mmseqs2 output.")

            if self.reverse_mmseqs2:
                col_list = [
                    "gene", "mmseqs_phrog", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal",
                    "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
                ]
            else:
                col_list = [
                    "mmseqs_phrog", "gene", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal",
                    "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
                ]

            mmseqs_df = pl.read_csv(
                mmseqs_file, separator="\t", has_header=False, new_columns=col_list, infer_schema=False
            )

            # get best hit per gene (lowest e-value)
            tophits_df = (
                mmseqs_df
                .with_columns(pl.col("mmseqs_eVal").cast(pl.Float64))
                .sort("mmseqs_eVal")
                .unique(subset=["gene"], keep="first", maintain_order=True)
                .select(["gene", "mmseqs_phrog", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal"])
            )

            # get all genes by left-joining with the all gene df
            tophits_df = (
                all_genes_df
                .join(tophits_df, on="gene", how="left")
                .with_columns([
                    pl.col("mmseqs_phrog").fill_null("No_MMseqs"),
                    pl.col("mmseqs_alnScore").fill_null("No_MMseqs"),
                    pl.col("mmseqs_seqIdentity").fill_null("No_MMseqs"),
                    pl.col("mmseqs_eVal").fill_null("No_MMseqs").cast(pl.Utf8),
                ])
            )
        else:
            tophits_df = all_genes_df.with_columns([
                pl.lit("No_MMseqs").alias("mmseqs_phrog"),
                pl.lit("No_MMseqs").alias("mmseqs_alnScore"),
                pl.lit("No_MMseqs").alias("mmseqs_seqIdentity"),
                pl.lit("No_MMseqs").alias("mmseqs_eVal"),
            ])

        ####################
        # combine phrogs
        ####################

        if self.hmm_flag is True:
            tophits_df = process_pyhmmer_results(tophits_df, self.pyhmmer_results_dict)
        else:
            tophits_df = tophits_df.with_columns([
                pl.lit("No_HMM").alias("pyhmmer_phrog"),
                pl.lit("No_HMM").alias("pyhmmer_bitscore"),
                pl.lit("No_HMM").alias("pyhmmer_evalue"),
            ])

        # strip off phrog_
        tophits_df = tophits_df.with_columns([
            pl.col("mmseqs_phrog").str.replace("phrog_", ""),
            pl.col("pyhmmer_phrog").str.replace("phrog_", ""),
        ])

        ############
        # code to create 1 overall phrog column
        if self.mmseqs_flag is True:
            tophits_df = tophits_df.with_columns(pl.col("mmseqs_phrog").alias("phrog"))
            # add pyhmmer phrog for any entry without mmseqs (null after left join)
            tophits_df = tophits_df.with_columns(
                pl.when(
                    pl.col("phrog").is_null() & (pl.col("pyhmmer_phrog") != "No_PHROG")
                ).then(pl.col("pyhmmer_phrog"))
                .otherwise(pl.col("phrog"))
                .alias("phrog")
            )
        else:
            tophits_df = tophits_df.with_columns(pl.col("pyhmmer_phrog").alias("phrog"))

        # read in phrog annotation file
        phrog_annot_df = pl.read_csv(
            os.path.join(self.db_dir, "phrog_annot_v4.tsv"), separator="\t", infer_schema=False
        )
        phrog_annot_df = phrog_annot_df.with_columns(pl.col("phrog").cast(pl.Utf8))

        # get only the contig id not the full description from mmseqs2 output
        tophits_df = tophits_df.with_columns(
            pl.col("gene").str.split(" ").list.first().alias("gene")
        )

        tophits_df = tophits_df.join(phrog_annot_df, on="phrog", how="left")
        tophits_df = tophits_df.with_columns([
            pl.col("annot").str.replace("No_PHROG", "hypothetical protein"),
            pl.col("category").str.replace("No_PHROG", "unknown function"),
        ])

        # cast mmseqs columns to str and fill nulls
        for col in ["mmseqs_phrog", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal", "color"]:
            tophits_df = tophits_df.with_columns(
                pl.col(col).cast(pl.Utf8).fill_null("No_PHROG")
            )

        # get phrog as string
        tophits_df = tophits_df.with_columns(pl.col("phrog").cast(pl.Utf8))

        # drop and re-merge for final annotation
        tophits_df = tophits_df.drop(["color", "annot", "category"])
        tophits_df = tophits_df.join(phrog_annot_df, on="phrog", how="left")
        tophits_df = tophits_df.with_columns([
            pl.col("annot").fill_null("hypothetical protein"),
            pl.col("category").fill_null("unknown function"),
            pl.col("color").fill_null("None"),
        ])
        # phrog_annot_v4.tsv stores "NA" as a literal string for unannotated phrogs.
        # pandas treated "NA" as NaN so fillna("hypothetical protein") worked;
        # polars infer_schema=False keeps "NA" as a string, so fill_null has no effect.
        # Replace explicitly here to match the old pandas behaviour.
        # Example: "NA" (old, pandas NaN-filled) → "hypothetical protein" (new, explicit).
        tophits_df = tophits_df.with_columns([
            pl.when(pl.col("annot") == "NA")
              .then(pl.lit("hypothetical protein"))
              .otherwise(pl.col("annot"))
              .alias("annot"),
            pl.when(pl.col("category") == "NA")
              .then(pl.lit("unknown function"))
              .otherwise(pl.col("category"))
              .alias("category"),
        ])

        ##### length df
        fasta_sequences = SeqIO.parse(open(self.input_fasta), "fasta")
        contig_names = []
        lengths = []
        for fasta in fasta_sequences:
            contig_names.append(fasta.id)
            lengths.append(len(fasta.seq))

        length_df = pl.DataFrame({"gene": contig_names, "length": lengths})
        self.length_df = length_df

        # left join length to tophits
        tophits_df = length_df.join(tophits_df, on="gene", how="left")

        # process vfdb results
        (tophits_df, vfdb_results) = process_vfdb_results(
            self.out_dir, tophits_df, proteins_flag=True, reverse_mmseqs2=self.reverse_mmseqs2
        )
        # process CARD results
        (tophits_df, card_results) = process_card_results(
            self.out_dir, tophits_df, self.db_dir, proteins_flag=True, reverse_mmseqs2=self.reverse_mmseqs2
        )

        # Rename the "gene" column to "ID"
        tophits_df = tophits_df.rename({"gene": "ID"})

        # Reorder columns: ID, length, then phrog, annot, category in positions 2,3,4
        cols = tophits_df.columns
        desired_order = [col for col in cols if col not in ["phrog", "annot", "category"]]
        desired_order.insert(2, "phrog")
        desired_order.insert(3, "annot")
        desired_order.insert(4, "category")
        tophits_df = tophits_df.select(desired_order)

        # fill any remaining nulls
        null_fills = {
            "phrog": "No_PHROG",
            "annot": "hypothetical protein",
            "category": "unknown function",
            "mmseqs_phrog": "No_MMseqs",
            "mmseqs_alnScore": "No_MMseqs",
            "mmseqs_seqIdentity": "No_MMseqs",
            "mmseqs_eVal": "No_MMseqs",
            "pyhmmer_phrog": "No_PHROGs_HMM",
            "pyhmmer_bitscore": "No_PHROGs_HMM",
            "pyhmmer_evalue": "No_PHROGs_HMM",
            "color": "None",
        }
        for col, fill_val in null_fills.items():
            if col in tophits_df.columns:
                tophits_df = tophits_df.with_columns(pl.col(col).fill_null(fill_val))

        tophits_df.write_csv(
            os.path.join(self.out_dir, f"{self.prefix}_full_merged_output.tsv"),
            separator="\t",
        )

        self.tophits_df = tophits_df
        self.vfdb_results = vfdb_results
        self.card_results = card_results

        summary_df = self.tophits_df.select(["ID", "length", "phrog", "annot", "category"])
        summary_df.write_csv(
            os.path.join(self.out_dir, f"{self.prefix}_summary_output.tsv"),
            separator="\t",
        )

    def update_fasta_headers(self):
        """
        updates FASTA header with description
        """
        fasta_output_aas = os.path.join(self.out_dir, f"{self.prefix}.faa")
        annots = self.tophits_df["annot"].to_list()

        with open(fasta_output_aas, "w") as aa_fa:
            for i, aa_record in enumerate(SeqIO.parse(self.input_fasta, "fasta")):
                aa_record.description = (
                    str(aa_record.description)
                    + " "
                    + str(annots[i])
                )
                SeqIO.write(aa_record, aa_fa, "fasta")


def main():
    """Entry point for `pharokka proteins` command."""
    # The logger.error → sys.exit(1) sink is registered once per process
    # by the entry-point dispatcher (cli.main or pharokka_scripts._legacy_shim).
    args = get_input_proteins()

    if args.citation is True:
        logger.info("If you use Pharokka in your research, please cite:")
        logger.info(
            "George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde."
        )
        logger.info("Pharokka: a fast scalable bacteriophage annotation tool.")
        logger.info("Bioinformatics, Volume 39, Issue 1, January 2023, btac776.")
        logger.info("https://doi.org/10.1093/bioinformatics/btac776.")
        logger.info(
            "You should also cite the full list of tools Pharokka uses, which can be found at https://github.com/gbouras13/pharokka#citation."
        )
        sys.exit()

    # set the prefix
    if args.prefix == "Default":
        prefix = "pharokka_proteins"
    else:
        prefix = args.prefix

    # instantiate outdir
    out_dir = instantiate_dirs(args.outdir, False, args.force)  # args.meta always false

    # set log dir
    logdir = Path(f"{out_dir}/logs")

    # set the db dir
    if args.database == "Default":
        _here = os.path.dirname(os.path.realpath(__file__))
        _src = os.path.dirname(_here)
        _root = os.path.dirname(_src)
        db_dir = os.path.join(_root, "databases/")
    else:
        db_dir = args.database

    # get start time
    start_time = time.time()
    # initial logging stuff
    log_file = os.path.join(args.outdir, f"pharokka_proteins_{start_time}.log")
    # adds log file
    logger.add(log_file)

    # preamble
    logger.info(f"Starting Pharokka v{get_version()}")
    logger.info("Running pharokka proteins to annotate proteins.")
    logger.info("Command executed: {}", args)
    logger.info("Repository homepage is https://github.com/gbouras13/pharokka")
    logger.info("Written by George Bouras: george.bouras@adelaide.edu.au")

    # check the database is installed
    logger.info("Checking database installation.")
    database_installed = check_db_installation(db_dir)
    if database_installed is True:
        logger.info("All databases have been successfully checked.")
    else:
        logger.error(
            "The database directory was unsuccessfully checked. Please run pharokka install"
        )

    # dependencies
    logger.info("Checking dependencies.")
    check_dependencies(False)  # to check pharokka proteins, don't need mash

    # instantiation/checking fasta
    validate_fasta(args.infile)
    validate_threads(args.threads)

    # define input
    input_fasta = args.infile

    ########
    # mmseqs2 and hmm decisions
    ########

    # can't have hmm_only and mmseqs2_only
    if args.hmm_only is True and args.mmseqs2_only is True:
        logger.error(
            "You have specified --fast or --hmm_only and --mmseqs2_only. This is impossible. Please choose one or the other."
        )

    # by default
    mmseqs_flag = True
    hmm_flag = True

    if args.hmm_only is True:
        logger.info("You have specified --hmm_only. MMseqs2 will not be run.")
        mmseqs_flag = False
    if args.mmseqs2_only is True:
        logger.info("You have specified --mmseqs2_only. Pyhmmer will not be run.")
        hmm_flag = False

    # running mmseqs2 on the 3 databases
    if mmseqs_flag is True:
        logger.info("Starting MMseqs2.")
        run_mmseqs_proteins(
            input_fasta,
            db_dir,
            out_dir,
            args.threads,
            logdir,
            args.evalue,
            args.reverse_mmseqs2,
            args.sensitivity,
            db_name="PHROG",
        )
        run_mmseqs_proteins(
            input_fasta,
            db_dir,
            out_dir,
            args.threads,
            logdir,
            args.evalue,
            args.reverse_mmseqs2,
            args.sensitivity,
            db_name="CARD",
        )
        run_mmseqs_proteins(
            input_fasta,
            db_dir,
            out_dir,
            args.threads,
            logdir,
            args.evalue,
            args.reverse_mmseqs2,
            args.sensitivity,
            db_name="VFDB",
        )

    if hmm_flag is True:
        # runs pyhmmer on PHROGs
        logger.info("Running PyHMMER.")
        best_results_pyhmmer = run_pyhmmer_proteins(
            input_fasta, db_dir, args.threads, args.evalue
        )

    #################################################
    # post processing
    #################################################

    logger.info("Post Processing Output.")

    # instanatiate the class with some of the params
    pharok = Pharok_Prot()
    pharok.out_dir = out_dir
    pharok.db_dir = db_dir
    pharok.prefix = prefix
    pharok.input_fasta = input_fasta
    pharok.mmseqs_flag = mmseqs_flag
    pharok.reverse_mmseqs2 = args.reverse_mmseqs2
    pharok.hmm_flag = hmm_flag
    if pharok.hmm_flag is True:
        pharok.pyhmmer_results_dict = best_results_pyhmmer

    # post process results
    pharok.process_dataframes()

    # updates fasta headers
    pharok.update_fasta_headers()

    # cleanup
    remove_directory(os.path.join(out_dir, "mmseqs"))
    remove_directory(os.path.join(out_dir, "CARD"))
    remove_directory(os.path.join(out_dir, "vfdb"))
    remove_directory(os.path.join(out_dir, "VFDB_dir"))
    remove_directory(os.path.join(out_dir, "tmp_dir"))
    remove_file(os.path.join(out_dir, "mmseqs_results.tsv"))
    remove_file(os.path.join(out_dir, "CARD_results.tsv"))
    remove_file(os.path.join(out_dir, "vfdb_results.tsv"))

    # Determine elapsed time
    elapsed_time = time.time() - start_time
    elapsed_time = round(elapsed_time, 2)

    # Show elapsed time for the process
    logger.info("Pharokka has finished")
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")

    logger.info("If you use Pharokka in your research, please cite:")
    logger.info(
        "George Bouras, Roshan Nepal, Ghais Houtak, Alkis James Psaltis, Peter-John Wormald, Sarah Vreugde."
    )
    logger.info("Pharokka: a fast scalable bacteriophage annotation tool.")
    logger.info("Bioinformatics, Volume 39, Issue 1, January 2023, btac776.")
    logger.info("https://doi.org/10.1093/bioinformatics/btac776.")
    logger.info(
        "You should also cite the full list of tools Pharokka uses, which can be found at https://github.com/gbouras13/pharokka#citation."
    )


if __name__ == "__main__":
    main()
