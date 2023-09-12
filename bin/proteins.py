"""
Holds functions for pharokka_proteins.py
"""

import argparse
import collections
import os
from argparse import RawTextHelpFormatter
from cmath import nan
from pathlib import Path
from re import T

import numpy as np
import pandas as pd
import pyhmmer
from Bio import SeqIO
from Bio.SeqUtils import GC
from external_tools import ExternalTool
from loguru import logger
from post_processing import (process_card_results, process_pyhmmer_results,
                             process_vfdb_results)
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile
from util import (count_contigs, get_contig_headers, get_version,
                  remove_directory)

Result = collections.namedtuple("Result", ["protein", "phrog", "bitscore", "evalue"])


pd.options.mode.chained_assignment = None


def get_input_proteins():
    parser = argparse.ArgumentParser(
        description="pharokka_proteins.py: fast phage protein annotation script",
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


def run_mmseqs_proteins(input_fasta, db_dir, out_dir, threads, logdir, evalue, db_name):
    """
    Runs mmseqs2 for pharokka_proteins.py
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
        tmp_dir = os.path.join(out_dir, "VFDB_dir/")
        profile_db = os.path.join(db_dir, "CARD")
        mmseqs_result_tsv = os.path.join(out_dir, "CARD_results.tsv")

    target_seqs = os.path.join(target_db_dir, "target_seqs")

    # make dir for target db
    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    # creates db for input
    mmseqs_createdb = ExternalTool(
        tool="mmseqs createdb",
        input=f"",
        output=f"{target_seqs}",
        params=f"{input_fasta}",  # param goes before output and mmseqs2 required order
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


# from itertools import chain


def run_pyhmmer_proteins(input_fasta, db_dir, threads, evalue):
    """
    Runs phymmer on phrogs
    :param input_fasta: input FASTA
    :param db_dir: database path
    :params threads: threads
    :param evalue: evalue threshold for pyhmmer
    :return: best_results - dictionary of top HMM hits for each protein
    """

    # define result
    Result = collections.namedtuple(
        "Result", ["protein", "phrog", "bitscore", "evalue"]
    )

    # run hmmscan and get all results
    results = []
    with pyhmmer.plan7.HMMFile(os.path.join(db_dir, "all_phrogs.h3m")) as hmms:  # hmms
        with pyhmmer.easel.SequenceFile(
            input_fasta, digital=True
        ) as seqs:  # amino acid sequences
            for hits in pyhmmer.hmmer.hmmscan(
                seqs, hmms, cpus=int(threads), E=float(evalue)
            ):  # run hmmscan
                protein = hits.query_name.decode()  # get protein from the hit
                for hit in hits:
                    if hit.included:
                        # include the hit to the result collection
                        results.append(
                            Result(protein, hit.name.decode(), hit.score, hit.evalue)
                        )

    # get  best results for each protein
    best_results = {}
    keep_protein = set()
    for result in results:
        if result.protein in best_results:
            previous_bitscore = best_results[result.protein].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.protein] = result
                keep_protein.add(result.protein)
            elif result.bitscore == previous_bitscore:
                if best_results[result.protein].phrog != hit.phrog:
                    keep_protein.remove(result.protein)
        else:
            best_results[result.protein] = result
            keep_protein.add(result.protein)

    return best_results


class Pharok_Prot:
    """pharokka_proteins.py Output Class"""

    def __init__(
        self,
        out_dir: Path = "output/",
        db_dir: Path = "db_dir/",
        prefix: str = "pharokka",
        input_fasta: Path = "input.fasta",
        pyhmmer_results_dict: dict = {
            "p1": Result(protein="p1", phrog="phrog_1", bitscore=1, evalue=2.01e-01)
        },
        tophits_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        vfdb_results: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        card_results: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        length_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        mmseqs_flag: bool = True,
        hmm_flag: bool = True,
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
        tophits_df: pd.DataFrame, required
            merged dataframe output
        vfdb_results: pd.DataFrame,
            vfdb dataframe output
        card_results: pd.DataFrame,
            CARD dataframe output
        length_df: pd.DataFrame,
            dataframe with lengths for each input contig
        run_mmseqs: bool
            whether MMseqs2 was run
        run_hmm: bool
            whether HMM was run
        """
        self.out_dir = out_dir
        self.db_dir = db_dir
        self.prefix = prefix
        self.input_fasta = input_fasta
        self.pyhmmer_results_dict = pyhmmer_results_dict
        self.tophits_df = tophits_df
        self.vfdb_results = vfdb_results
        self.card_results = card_results
        self.length_df = length_df
        self.mmseqs_flag = mmseqs_flag
        self.hmm_flag = hmm_flag

    def process_dataframes(self):
        """
        Processes and combines PHROGS, CARD and VFDB mmseqs output
        :param self: Pharok class object
        :param db_dir: database directory path
        :param out_dir: output directory path
        :return:
        """

        ####### ####### ####### ####### #######
        ####### ####### ####### ####### #######
        # get tophits
        ####### ####### ####### ####### #######
        ####### ####### ####### ####### #######
        if self.mmseqs_flag is True:
            # MMseqs PHROGs file
            mmseqs_file = os.path.join(self.out_dir, "mmseqs_results.tsv")
            logger.info("Processing mmseqs2 output.")
            col_list = [
                "mmseqs_phrog",
                "gene",
                "mmseqs_alnScore",
                "mmseqs_seqIdentity",
                "mmseqs_eVal",
                "qStart",
                "qEnd",
                "qLen",
                "tStart",
                "tEnd",
                "tLen",
            ]
            mmseqs_df = pd.read_csv(
                mmseqs_file, delimiter="\t", index_col=False, names=col_list
            )
            # get list of genes
            genes = mmseqs_df.gene.unique()

            # instantiate tophits list
            tophits = []

            for gene in genes:
                tmp_df = (
                    mmseqs_df.loc[mmseqs_df["gene"] == gene]
                    .sort_values("mmseqs_eVal")
                    .reset_index(drop=True)
                    .loc[0]
                )
                tophits.append(
                    [
                        tmp_df.gene,
                        tmp_df.mmseqs_phrog,
                        tmp_df.mmseqs_alnScore,
                        tmp_df.mmseqs_seqIdentity,
                        tmp_df.mmseqs_eVal,
                    ]
                )

            # create tophits df
            tophits_df = pd.DataFrame(
                tophits,
                columns=[
                    "gene",
                    "mmseqs_phrog",
                    "mmseqs_alnScore",
                    "mmseqs_seqIdentity",
                    "mmseqs_eVal",
                ],
            )

        # tophits for pyhmmer only
        else:
            prot_count = count_contigs(self.input_fasta)
            ids = get_contig_headers(self.input_fasta)
            # create tophits df
            # Create a dictionary with the column names and their corresponding values
            data = {
                "gene": ids,
                "mmseqs_phrog": ["No_MMseqs"] * prot_count,
                "mmseqs_alnScore": ["No_MMseqs"] * prot_count,
                "mmseqs_seqIdentity": ["No_MMseqs"] * prot_count,
                "mmseqs_eVal": ["No_MMseqs"] * prot_count,
            }

            # Create a DataFrame from the dictionary
            tophits_df = pd.DataFrame(data)

        if len(tophits_df["mmseqs_phrog"]) == 0:
            tophits_df["mmseqs_top_hit"] = "No_PHROG"
        else:
            if self.mmseqs_flag is True:  # trim the rubbish if mmseqs2 is on
                tophits_df[["mmseqs_phrog", "mmseqs_top_hit"]] = tophits_df[
                    "mmseqs_phrog"
                ].str.split(" ## ", expand=True)
            else:  # no mmseqs2 hits
                tophits_df["mmseqs_top_hit"] = "No_MMseqs_PHROG_hit"

        ####################
        # combine phrogs
        ####################

        # Adds pyhmmer results if tru
        if self.hmm_flag is True:
            tophits_df = process_pyhmmer_results(tophits_df, self.pyhmmer_results_dict)
        else:
            tophits_df["pyhmmer_phrog"] = "No_HMM"
            tophits_df["pyhmmer_bitscore"] = "No_HMM"
            tophits_df["pyhmmer_evalue"] = "No_HMM"

        # strip off phrog_ for both
        tophits_df["mmseqs_phrog"] = tophits_df["mmseqs_phrog"].str.replace(
            "phrog_", ""
        )
        tophits_df["pyhmmer_phrog"] = tophits_df["pyhmmer_phrog"].str.replace(
            "phrog_", ""
        )

        ############
        # code to create 1 overall phrog column
        # pick the mmseqs PHROG column first if it was run

        if self.mmseqs_flag is True:
            tophits_df["phrog"] = tophits_df["mmseqs_phrog"]
            # add pyhmmer phrog for any entry without mmseqs
            for index, row in tophits_df.iterrows():
                if isinstance(
                    row["phrog"], float
                ):  # for all the rows without an mmseqs2 PHROG will be nan - floats
                    # to write all hits where there was no mmseqs but there was a hmm
                    if row["pyhmmer_phrog"] != "No_PHROG":
                        tophits_df.at[index, "phrog"] = row["pyhmmer_phrog"]

        else:  # only need to worry about pyhmmer
            tophits_df["phrog"] = tophits_df["pyhmmer_phrog"]

        # read in phrog annotaion file
        phrog_annot_df = pd.read_csv(
            os.path.join(self.db_dir, "phrog_annot_v4.tsv"), sep="\t", index_col=False
        )
        phrog_annot_df["phrog"] = phrog_annot_df["phrog"].astype(str)

        # merge phrog
        tophits_df = tophits_df.merge(phrog_annot_df, on="phrog", how="left")
        tophits_df = tophits_df.replace(np.nan, "No_PHROG", regex=True)
        # convert no phrog to hyp protein
        tophits_df["annot"] = tophits_df["annot"].str.replace(
            "No_PHROG", "hypothetical protein"
        )
        tophits_df["category"] = tophits_df["category"].str.replace(
            "No_PHROG", "unknown function"
        )

        # # replace with No_PHROG if nothing found
        tophits_df.loc[
            tophits_df["mmseqs_phrog"] == "No_PHROG", "mmseqs_phrog"
        ] = "No_PHROG"
        tophits_df.loc[
            tophits_df["mmseqs_alnScore"] == "No_PHROG", "mmseqs_alnScore"
        ] = "No_PHROG"
        tophits_df.loc[
            tophits_df["mmseqs_seqIdentity"] == "No_PHROG", "mmseqs_seqIdentity"
        ] = "No_PHROG"
        tophits_df.loc[
            tophits_df["mmseqs_eVal"] == "No_PHROG", "mmseqs_eVal"
        ] = "No_PHROG"
        tophits_df.loc[
            tophits_df["mmseqs_top_hit"] == "No_PHROG", "mmseqs_top_hit"
        ] = "No_PHROG"
        tophits_df.loc[tophits_df["color"] == "No_PHROG", "color"] = "No_PHROG"

        # get phrog
        tophits_df["phrog"] = tophits_df["phrog"].astype(str)

        # drop existing color annot category cols and
        tophits_df = tophits_df.drop(columns=["color", "annot", "category"])
        tophits_df = tophits_df.merge(phrog_annot_df, on="phrog", how="left")
        tophits_df["annot"] = tophits_df["annot"].replace(
            nan, "hypothetical protein", regex=True
        )
        tophits_df["category"] = tophits_df["category"].replace(
            nan, "unknown function", regex=True
        )
        tophits_df["color"] = tophits_df["color"].replace(nan, "None", regex=True)

        ##### length df

        fasta_sequences = SeqIO.parse(open(self.input_fasta), "fasta")
        contig_names = []
        lengths = []

        for fasta in fasta_sequences:
            contig_names.append(fasta.id)
            lengths.append(len(fasta.seq))

        length_df = pd.DataFrame(
            {
                "gene": contig_names,
                "length": lengths,
            }
        )
        self.length_df = length_df

        # merge the length df into the tophits
        tophits_df = length_df.merge(tophits_df, on="gene", how="left")

        # process vfdb results
        # handles empty files without a problem
        (tophits_df, vfdb_results) = process_vfdb_results(self.out_dir, tophits_df)
        # process CARD results
        (tophits_df, card_results) = process_card_results(
            self.out_dir, tophits_df, self.db_dir
        )

        # Rename the "gene" column to "id"
        tophits_df.rename(columns={"gene": "ID"}, inplace=True)

        desired_order = [
            col
            for col in tophits_df.columns
            if col not in ["phrog", "annot", "category"]
        ]
        desired_order.insert(2, "phrog")
        desired_order.insert(3, "annot")
        desired_order.insert(4, "category")

        # Creating a new DataFrame with columns rearranged
        tophits_df = tophits_df[desired_order]

        # save

        tophits_df.to_csv(
            os.path.join(self.out_dir, f"{self.prefix}_full_merged_output.tsv"),
            sep="\t",
            index=False,
        )

        # save results
        self.tophits_df = tophits_df
        self.vfdb_results = vfdb_results
        self.card_results = card_results

        summary_df = self.tophits_df[["ID", "length", "phrog", "annot", "category"]]

        # .tsv file with all data full_merged_output
        summary_df.to_csv(
            os.path.join(self.out_dir, f"{self.prefix}_summary_output.tsv"),
            sep="\t",
            index=False,
        )

    def update_fasta_headers(self):
        """
        updates FASTA header with description
        """
        self.input_fasta

        # define outputs
        fasta_output_aas = os.path.join(self.out_dir, f"{self.prefix}.faa")

        # amino acids

        with open(fasta_output_aas, "w") as aa_fa:
            i = 0
            for aa_record in SeqIO.parse(self.input_fasta, "fasta"):
                aa_record.description = (
                    str(aa_record.description)
                    + " "
                    + str(self.tophits_df["annot"].iloc[i])
                )
                SeqIO.write(aa_record, aa_fa, "fasta")
                i += 1
