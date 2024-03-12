import collections
import os
import random
import string
from cmath import nan
from pathlib import Path
from re import T

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from loguru import logger
from processes import convert_gff_to_gbk
from util import remove_directory, remove_file, touch_file

pd.options.mode.chained_assignment = None

Result = collections.namedtuple("Result", ["protein", "phrog", "bitscore", "evalue"])


class Pharok:
    """Pharokka Output Class"""

    def __init__(
        self,
        out_dir: Path = "output/",
        db_dir: Path = "db_dir/",
        prefix: str = "pharokka",
        gene_predictor: str = "phanotate",
        input_fasta: Path = "input.fasta",
        meta_mode: bool = False,
        locustag: str = "Random",
        pyhmmer_results_dict: dict = {
            "p1": Result(protein="p1", phrog="phrog_1", bitscore=1, evalue=2.01e-01)
        },
        custom_pyhmmer_results_dict: dict = {
            "p1": Result(protein="p1", phrog="phrog_1", bitscore=1, evalue=2.01e-01)
        },
        merged_df: pd.DataFrame() = pd.DataFrame(
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
        gff_df: pd.DataFrame() = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}),
        locus_df: pd.DataFrame() = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}),
        prot_seq_df: pd.DataFrame() = pd.DataFrame(
            {"col1": [1, 2, 3], "col2": [4, 5, 6]}
        ),
        tmrna_flag: bool = False,
        trna_empty: bool = False,
        crispr_count: int = 0,
        coding_table: int = 11,
        mmseqs_flag: bool = True,
        hmm_flag: bool = True,
        custom_hmm_flag: bool = False,
        phanotate_version: str = "1.5.0",
        pyrodigal_version: str = "3.0.0",
        pyrodigal_gv_version: str = "0.1.0",
        trna_version: str = "2.0.12",
        aragorn_version: str = "1.2.41",
        minced_version: str = "0.4.2",
        skip_extra_annotations: bool = False,
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
        gene_predictor: str, gene predictor used
            phanotate or prodigal
        input_fasta: Path
            input FASTA file with phage contigs
        meta_mode: bool
            boolean whether pharokka is being run in meta mode
        locustag: str
            locustag. Will be 'Random' if the user has not specified a locus tag
        pyhmmer_results_dict: dict
            dictionary of Result tuples from pyhmmer.
        merged_df: pd.DataFrame, required
            merged dataframe output
        vfdb_results: pd.DataFrame,
            vfdb dataframe output
        card_results: pd.DataFrame,
            CARD dataframe output
        length_df: pd.DataFrame,
            dataframe with lengths for each input contig
        gff_df: pd.DataFrame,
            dataframe with the output for the gff
        locus_df: pd.DataFrame,
            dataframe with locus tag information
        tmrna_flag: boolean,
            flag denoting whether there is a tmRNA in the output (per aragorn)
        trna_empty: boolean,
            flag denoting whether the trna file is empty (True means no trnas)
        crispr_count: int,
            number of CRISPRs
        coding_table: int.
            number denoting the prodigal coding table (default 11)
        mmseqs_flag: bool
            whether MMseqs2 was run
        hmm_flag: bool
            whether HMM was run
        custom_hmm_flag: bool
            whether a custom db of HMMs was run
        phanotate_version: str
            phanotate_version from check_dependencies()
        prodigal_version: str
            prodigal_version from check_dependencies()
        pyrodigal_gv_version: str
            pyrodigal_gv_version from check dependencies()
        trna_version: str
            trnascan_version from check_dependencies()
        aragorn_version: str
            aragorn_version from check_dependencies()
        minced_version: str
            minced_version from check dependencies()
        prot_seq_df: pd.DataFrame,
            dataframe with protein sequence  information for each egene
        skip_extra_annotations: bool
            boolean whether extra annotations are skipped
        """
        self.out_dir = out_dir
        self.db_dir = db_dir
        self.prefix = prefix
        self.gene_predictor = gene_predictor
        self.input_fasta = input_fasta
        self.meta_mode = meta_mode
        self.locustag = locustag
        self.pyhmmer_results_dict = pyhmmer_results_dict
        self.custom_pyhmmer_results_dict = custom_pyhmmer_results_dict
        self.merged_df = merged_df
        self.vfdb_results = vfdb_results
        self.card_results = card_results
        self.length_df = length_df
        self.gff_df = gff_df
        self.locus_df = locus_df
        self.tmrna_flag = tmrna_flag
        self.trna_empty = trna_empty
        self.crispr_count = crispr_count
        self.coding_table = coding_table
        self.mmseqs_flag = mmseqs_flag
        self.hmm_flag = hmm_flag
        self.custom_hmm_flag = custom_hmm_flag
        self.phanotate_version = phanotate_version
        self.pyrodigal_version = pyrodigal_version
        self.pyrodigal_gv_version = pyrodigal_gv_version
        self.trna_version = trna_version
        self.aragorn_version = aragorn_version
        self.minced_version = minced_version
        self.prot_seq_df = prot_seq_df
        self.skip_extra_annotations = skip_extra_annotations

    def process_results(self):
        """
        Processes and combines PHROGS, CARD and VFDB mmseqs output
        :param self: Pharok class object
        :param db_dir: database directory path
        :param out_dir: output directory path
        :param gene_predictor: CDS predictor (phanotate or prodigal)
        :return:
        """

        # left join mmseqs top hits to cds df
        # read in the cds cdf
        cds_file = os.path.join(self.out_dir, "cleaned_" + self.gene_predictor + ".tsv")

        col_list = ["start", "stop", "frame", "contig", "score", "gene"]
        dtype_dict = {
            "start": int,
            "stop": int,
            "frame": str,
            "contig": str,
            "score": str,
            "gene": str,
        }

        cds_df = pd.read_csv(
            cds_file,
            sep="\t",
            index_col=False,
            names=col_list,
            dtype=dtype_dict,
            skiprows=1,
        )
        cds_df["contig"] = cds_df["contig"].astype(str)

        ###########################################
        # add the sequence to the df for the genbank conversion later on
        fasta_input_aas_tmp = os.path.join(
            self.out_dir, f"{self.gene_predictor}_aas_tmp.fasta"
        )
        prot_dict = SeqIO.to_dict(SeqIO.parse(fasta_input_aas_tmp, "fasta"))

        # make a copy of cds_df
        self.prot_seq_df = cds_df.copy()

        # to match the output for gff
        self.prot_seq_df[["gene", "st"]] = self.prot_seq_df["gene"].str.split(
            " ", expand=True
        )

        self.prot_seq_df = self.prot_seq_df.drop(columns=["st"])

        # get sequences for each gene in df
        self.prot_seq_df["sequence"] = "MA"

        for index, row in self.prot_seq_df.iterrows():
            # get the gene id
            gene = row["gene"]
            if gene in prot_dict.keys():
                # add the AA sequence
                self.prot_seq_df.at[index, "sequence"] = prot_dict[gene].seq

        ##########################################
        # create the tophits_df and write it to file
        if self.mmseqs_flag is True:
            tophits_df = create_mmseqs_tophits(self.out_dir)

        else:
            # create tophits df
            # Create a dictionary with the column names and their corresponding values
            data = {
                "mmseqs_phrog": ["No_MMseqs"] * len(cds_df),
                "gene": cds_df["gene"],
                "mmseqs_alnScore": ["No_MMseqs"] * len(cds_df),
                "mmseqs_seqIdentity": ["No_MMseqs"] * len(cds_df),
                "mmseqs_eVal": ["No_MMseqs"] * len(cds_df),
            }

            # Create a DataFrame from the dictionary
            tophits_df = pd.DataFrame(data)

        # convert the gene to string for the merge
        cds_df["gene"] = cds_df["gene"].astype(str)
        tophits_df["gene"] = tophits_df["gene"].astype(str)
        cds_df = cds_df[cds_df["start"].notna()]
        cds_df = cds_df.dropna()

        # merge top hits into the cds df
        merged_df = cds_df.merge(tophits_df, on="gene", how="left")

        # get best protein from top mmseqs2 hit
        # add test if empty - crashes if no gene call hits

        if len(tophits_df["mmseqs_phrog"]) == 0:
            merged_df["mmseqs_top_hit"] = "No_PHROG"
        else:
            if self.mmseqs_flag is True:  # trim the rubbish if mmseqs2 is on
                merged_df[["mmseqs_phrog", "mmseqs_top_hit"]] = merged_df[
                    "mmseqs_phrog"
                ].str.split(" ## ", expand=True)
            else:  # no mmseqs2 hits
                merged_df["mmseqs_top_hit"] = "No_MMseqs_PHROG_hit"

        ####################
        # combine phrogs
        ####################

        # Adds pyhmmer results if tru
        if self.hmm_flag is True:
            merged_df = process_pyhmmer_results(merged_df, self.pyhmmer_results_dict)
        else:
            merged_df["pyhmmer_phrog"] = "No_PHROGs_HMM"
            merged_df["pyhmmer_bitscore"] = "No_PHROGs_HMM"
            merged_df["pyhmmer_evalue"] = "No_PHROGs_HMM"

        # add custom db results
        if self.custom_hmm_flag is True:
            merged_df = process_custom_pyhmmer_results(
                merged_df, self.custom_pyhmmer_results_dict
            )
        else:
            merged_df["custom_hmm_id"] = "No_custom_HMM"
            merged_df["custom_hmm_bitscore"] = "No_custom_HMM"
            merged_df["custom_hmm_evalue"] = "No_custom_HMM"

        # strip off phrog_ for both
        merged_df["mmseqs_phrog"] = merged_df["mmseqs_phrog"].str.replace("phrog_", "")
        merged_df["pyhmmer_phrog"] = merged_df["pyhmmer_phrog"].str.replace(
            "phrog_", ""
        )

        ############
        # code to create 1 overall phrog column
        # pick the mmseqs PHROG column first if it was run

        if self.mmseqs_flag is True:
            merged_df["phrog"] = merged_df["mmseqs_phrog"]
            # add pyhmmer phrog for any entry without mmseqs
            for index, row in merged_df.iterrows():
                if isinstance(
                    row["phrog"], float
                ):  # for all the rows without an mmseqs2 PHROG will be nan - floats
                    # to write all hits where there was no mmseqs but there was a hmm
                    if row["pyhmmer_phrog"] != "No_PHROGs_HMM":
                        merged_df.at[index, "phrog"] = row["pyhmmer_phrog"]

        else:  # only need to worry about pyhmmer
            merged_df["phrog"] = merged_df["pyhmmer_phrog"]

        # read in phrog annotaion file
        phrog_annot_df = pd.read_csv(
            os.path.join(self.db_dir, "phrog_annot_v4.tsv"), sep="\t", index_col=False
        )
        phrog_annot_df["phrog"] = phrog_annot_df["phrog"].astype(str)

        # merge phrog
        merged_df = merged_df.merge(phrog_annot_df, on="phrog", how="left")
        merged_df = merged_df.replace(np.nan, "No_PHROG", regex=True)
        # convert no phrog to hyp protein
        merged_df["annot"] = merged_df["annot"].str.replace(
            "No_PHROG", "hypothetical protein"
        )
        merged_df["category"] = merged_df["category"].str.replace(
            "No_PHROG", "unknown function"
        )

        # add columns
        if self.gene_predictor == "phanotate":
            merged_df["Method"] = f"PHANOTATE_{self.phanotate_version}"
        elif self.gene_predictor == "prodigal":
            merged_df["Method"] = f"Pyrodigal_{self.pyrodigal_version}"
        elif self.gene_predictor == "genbank":
            merged_df["Method"] = "CUSTOM"
        elif self.gene_predictor == "prodigal-gv":
            merged_df["Method"] = f"Pyrodigal-gv_{self.pyrodigal_gv_version}"
        merged_df["Region"] = "CDS"

        # # replace with No_PHROG if nothing found
        merged_df.loc[
            merged_df["mmseqs_phrog"] == "No_PHROG", "mmseqs_phrog"
        ] = "No_PHROG"
        merged_df.loc[
            merged_df["mmseqs_alnScore"] == "No_PHROG", "mmseqs_alnScore"
        ] = "No_PHROG"
        merged_df.loc[
            merged_df["mmseqs_seqIdentity"] == "No_PHROG", "mmseqs_seqIdentity"
        ] = "No_PHROG"
        merged_df.loc[
            merged_df["mmseqs_eVal"] == "No_PHROG", "mmseqs_eVal"
        ] = "No_PHROG"
        merged_df.loc[
            merged_df["mmseqs_top_hit"] == "No_PHROG", "mmseqs_top_hit"
        ] = "No_PHROG"
        merged_df.loc[merged_df["color"] == "No_PHROG", "color"] = "No_PHROG"

        # get phrog
        merged_df["phrog"] = merged_df["phrog"].astype(str)

        # drop existing color annot category cols and
        merged_df = merged_df.drop(columns=["color", "annot", "category"])
        merged_df = merged_df.merge(phrog_annot_df, on="phrog", how="left")
        merged_df["annot"] = merged_df["annot"].replace(
            nan, "hypothetical protein", regex=True
        )
        merged_df["category"] = merged_df["category"].replace(
            nan, "unknown function", regex=True
        )
        merged_df["color"] = merged_df["color"].replace(nan, "None", regex=True)

        # process vfdb results
        # handles empty files without a problem
        (merged_df, vfdb_results) = process_vfdb_results(
            self.out_dir, merged_df, proteins_flag=False
        )
        # process CARD results
        (merged_df, card_results) = process_card_results(
            self.out_dir, merged_df, self.db_dir, proteins_flag=False
        )

        self.merged_df = merged_df
        self.vfdb_results = vfdb_results
        self.card_results = card_results

    def get_contig_name_lengths(self):
        """
        Gets contig name and length in the input fasta file and calculates gc.
        Also adds translation table
        :param fasta_input: input fasta file
        :return: length_df a pandas dataframe (to class)
        """

        fasta_sequences = SeqIO.parse(open(self.input_fasta), "fasta")

        if self.gene_predictor == "prodigal-gv":
            # define col list
            col_list = [
                "contig",
                "Method",
                "Region",
                "start",
                "stop",
                "score",
                "frame",
                "phase",
                "attributes",
            ]
            # read gff (no fasta output)
            pyrodigal_gv_gff = pd.read_csv(
                os.path.join(self.out_dir, "prodigal-gv_out.gff"),
                delimiter="\t",
                index_col=False,
                names=col_list,
            )

            pyrodigal_gv_gff[["attributes", "transl_table"]] = pyrodigal_gv_gff[
                "attributes"
            ].str.split("transl_table=", expand=True)
            pyrodigal_gv_gff[["transl_table", "rest"]] = pyrodigal_gv_gff[
                "transl_table"
            ].str.split(";conf", expand=True)
            # drop and then remove duplicates in df
            pyrodigal_gv_gff = pyrodigal_gv_gff.drop(
                columns=[
                    "rest",
                    "Method",
                    "Region",
                    "start",
                    "stop",
                    "score",
                    "frame",
                    "phase",
                    "attributes",
                ]
            )
            # Remove duplicate rows based on all columns
            pyrodigal_gv_gff = pyrodigal_gv_gff.drop_duplicates()
            # Convert to a dictionary
            transl_table_dict = pyrodigal_gv_gff.set_index("contig")[
                "transl_table"
            ].to_dict()

        contig_names = []
        lengths = []
        gc = []
        transl_tables = []

        transl_table = "11"
        if self.gene_predictor == "phanotate":
            transl_table = "11"
        elif self.gene_predictor == "prodigal":
            transl_table = self.coding_table
        elif self.gene_predictor == "genbank":
            transl_table = "custom_gene_calls_from_genbank"

        for record in fasta_sequences:
            contig_names.append(record.id)
            lengths.append(len(record.seq))
            gc.append(round(gc_fraction(record.seq), 2))
            # pyrodigal-gv lookup from the dict
            if self.gene_predictor == "prodigal-gv":
                # try catch clause if contig too small to have a gene
                try:
                    transl_table = transl_table_dict[record.id]
                except:
                    transl_table = "No_CDS_called"

            transl_tables.append(transl_table)

        length_df = pd.DataFrame(
            {
                "contig": contig_names,
                "length": lengths,
                "gc_perc": gc,
                "transl_table": transl_tables,
            }
        )
        self.length_df = length_df

    def parse_aragorn(self):
        """
        Parses aragorn output file
        :param out_dir: output directory path
        :param prefix: prefix
        :param length_df: a pandas df as output from get_contig_name_lengths()
        :return:- to class - tmrna_flag Boolean whethere there is tmrna or not
        """
        aragorn_file = os.path.join(self.out_dir, self.prefix + "_aragorn.txt")
        f = open(aragorn_file)
        lines = f.readlines()
        contig_count = len(self.length_df["contig"])
        tmrna_flag = False
        contig_names = []
        methods = []
        regions = []
        starts = []
        stops = []
        scores = []
        frames = []
        phases = []
        attributes = []
        # if there is only one contig
        if contig_count == 1:
            # if no trnas
            if int(lines[1][0]) == 0:
                tmrna_df = pd.DataFrame(
                    {
                        "contig": "",
                        "Method": "",
                        "Region": "",
                        "start": "",
                        "stop": "",
                        "score": "",
                        "frame": "",
                        "phase": "",
                        "attributes": "",
                    },
                    index=[0],
                )
            else:
                tmrna_flag = True
                # get all lines with tmrnas
                tmrna_lines = lines[2:]
                for line in tmrna_lines:
                    split = line.split()
                    start_stops = split[2].replace("[", "").replace("]", "").split(",")
                    contig = self.length_df["contig"][0]
                    method = f"Aragorn_{self.aragorn_version}"
                    region = "tmRNA"
                    start = start_stops[0].replace(
                        "c", ""
                    )  # tmrna output is [start,stop] or c[start, stop] so need to remove c also for some phages
                    stop = start_stops[1]
                    score = "."
                    frame = "."
                    phase = "."
                    tag_peptide = split[3]
                    tag_peptide_seq = split[4]
                    attribute = (
                        "product=transfer-messenger RNA SsrA;tag_peptide="
                        + tag_peptide
                        + ";tag_peptide_sequence="
                        + tag_peptide_seq
                    )
                    contig_names.append(contig)
                    methods.append(method)
                    regions.append(region)
                    starts.append(start)
                    stops.append(stop)
                    scores.append(score)
                    frames.append(frame)
                    phases.append(phase)
                    attributes.append(attribute)
                tmrna_df = pd.DataFrame(
                    {
                        "contig": contig_names,
                        "Method": methods,
                        "Region": regions,
                        "start": starts,
                        "stop": stops,
                        "score": scores,
                        "frame": frames,
                        "phase": phases,
                        "attributes": attributes,
                    }
                )
        # two or more contigs
        else:
            i = 0  # line counter
            j = 0  # contig counter
            for line in lines:
                # contig meets these reqs
                if (
                    line[0] == ">" and "sensitivity" not in line
                ):  # sensititvity in the line at the end of the file
                    if "0 genes found" not in lines[i + 1]:
                        tmrna_flag = True
                        # number of trnas for this contig
                        tmrna_count = int(lines[i + 1][0])
                        # iterate over them
                        for k in range(tmrna_count):
                            tmrna_line = lines[i + 2 + k]
                            split = tmrna_line.split()
                            start_stops = (
                                split[2].replace("[", "").replace("]", "").split(",")
                            )
                            contig = self.length_df["contig"][j]
                            method = f"Aragorn_{self.aragorn_version}"
                            region = "tmRNA"
                            start = start_stops[0].replace(
                                "c", ""
                            )  # tmrna output is [start,stop] or c[start, stop] so need to remove c also for some phages
                            stop = start_stops[1]
                            score = "."
                            frame = "."
                            phase = "."
                            tag_peptide = split[3]
                            tag_peptide_seq = split[4]
                            attribute = (
                                "product=transfer-messenger RNA SsrA;tag_peptide="
                                + tag_peptide
                                + ";tag_peptide_sequence="
                                + tag_peptide_seq
                            )
                            contig_names.append(contig)
                            methods.append(method)
                            regions.append(region)
                            starts.append(start)
                            stops.append(stop)
                            scores.append(score)
                            frames.append(frame)
                            phases.append(phase)
                            attributes.append(attribute)
                    j += 1  # iterate contig
                # iterate line
                i += 1
            # write out the df
            tmrna_df = pd.DataFrame(
                {
                    "contig": contig_names,
                    "Method": methods,
                    "Region": regions,
                    "start": starts,
                    "stop": stops,
                    "score": scores,
                    "frame": frames,
                    "phase": phases,
                    "attributes": attributes,
                }
            )
        tmrna_df.to_csv(
            os.path.join(self.out_dir, self.prefix + "_aragorn.gff"),
            sep="\t",
            index=False,
            header=False,
        )
        self.tmrna_flag = tmrna_flag

    def create_gff(self):
        """
        Creates the pharokka.gff file
        :param merged_df: a pandas df as output from process_results() with results for cds
        :param length_df: a pandas df as output from get_contig_name_lengths()
        :param fasta_input: input fasta file
        :param out_dir: output directory path
        :param prefix: output prefix
        :param locustag: str whether or not to create a random locustag - will be 'Random' is so. Otherwise it is parsed
        :tmrna_flag boolean whether there are tmRNAs or not
        :return: locustag for the creation of the .tbl file, locus_df as df with consistent locus_tags


        merged_df,
        length_df,
        fasta_input,
        out_dir,
        prefix,
        locustag,
        tmrna_flag,
        meta_mode,

        """

        # create locus tag
        # locustag creation
        if self.locustag == "Random":
            # locus tag header 8 random letters
            self.locustag = "".join(
                random.choice(string.ascii_uppercase) for _ in range(8)
            )

        # get all contigs
        contigs = self.length_df["contig"].astype(str)
        self.merged_df["contig"] = self.merged_df["contig"].astype(str)

        # add the translation table
        # typing
        transl_table_df = self.length_df.drop(columns=["length", "gc_perc"])
        transl_table_df["contig"] = transl_table_df["contig"].astype(str)
        self.merged_df = self.merged_df.merge(transl_table_df, how="left", on="contig")

        ############ locus tag #########
        # write df for locus tag parsing
        # zfill - makes the CDS 4 digits trailing zeroes for vcontact
        # in meta mode,
        locus_df = self.merged_df
        subset_dfs = []

        # if meta mode is true
        if self.meta_mode == True:
            for contig in contigs:
                subset_df = locus_df[locus_df["contig"] == contig].reset_index()
                subset_df["count"] = subset_df.index
                # so not 0 indexed
                subset_df["count"] = subset_df["count"] + 1
                # z fill to make the locus tag 4
                subset_df["count"] = subset_df["count"].astype(str).str.zfill(4)
                subset_dfs.append(subset_df)
            locus_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
        else:
            locus_df["count"] = locus_df.index
            # so not 0 indexed
            locus_df["count"] = locus_df["count"] + 1
            # z fill to make the locus tag 4
            locus_df["count"] = locus_df["count"].astype(str).str.zfill(4)

        # get the locus tag
        if self.meta_mode == False:
            locus_df["locus_tag"] = self.locustag + "_CDS_" + locus_df["count"]
        else:  # for meta use contig names
            locus_df["locus_tag"] = locus_df.contig + "_CDS_" + locus_df["count"]

        # assign count and locus_tag to merged_df (for meta)
        self.merged_df["locus_tag"] = locus_df["locus_tag"]
        self.merged_df["count"] = locus_df["count"]

        #################################

        #########
        # rearrange start and stop so that start is always less than stop for gff
        #########

        cols = ["start", "stop"]
        # indices where start is greater than stop
        ixs = self.merged_df["frame"] == "-"
        # Where ixs is True, values are swapped
        self.merged_df.loc[ixs, cols] = (
            self.merged_df.loc[ixs, cols].reindex(columns=cols[::-1]).values
        )

        # set phase to be 0
        self.merged_df["phase"] = 0
        # create attributes
        # no custom

        self.merged_df["attributes"] = (
            "ID="
            + locus_df["locus_tag"].astype(str)
            + ";"
            + "transl_table="
            + locus_df["transl_table"].astype(str)
            + ";"
            + "phrog="
            + self.merged_df["phrog"].astype(str)
            + ";"
            + "top_hit="
            + self.merged_df["mmseqs_top_hit"].astype(str)
            + ";"
            + "locus_tag="
            + locus_df["locus_tag"].astype(str)
            + ";"
            + "function="
            + self.merged_df["category"].astype(str)
            + ";"
            + "product="
            + self.merged_df["annot"].astype(str)
        )
        # adds custom hmm database annotations
        self.merged_df.loc[
            self.merged_df["custom_hmm_id"] != "No_custom_HMM", "attributes"
        ] = (
            self.merged_df["attributes"].astype(str)
            + ";"
            + "custom_annotation="
            + self.merged_df["custom_hmm_id"].astype(str)
        )
        # adds VFDB
        self.merged_df.loc[
            self.merged_df["vfdb_short_name"] != "None", "attributes"
        ] = (
            self.merged_df["attributes"].astype(str)
            + ";"
            + "vfdb_short_name="
            + self.merged_df["vfdb_short_name"].astype(str)
            + ";"
            + "vfdb_description="
            + self.merged_df["vfdb_description"].astype(str)
            + ";"
            + "vfdb_species="
            + self.merged_df["vfdb_species"].astype(str)
        )
        # adds CARD
        self.merged_df.loc[
            self.merged_df["CARD_short_name"] != "None", "attributes"
        ] = (
            self.merged_df["attributes"].astype(str)
            + ";"
            + "CARD_short_name="
            + self.merged_df["CARD_short_name"].astype(str)
            + ";"
            + "AMR_Gene_Family="
            + self.merged_df["AMR_Gene_Family"].astype(str)
            + ";"
            + "CARD_species="
            + self.merged_df["CARD_species"].astype(str)
        )

        # save back
        # get gff dataframe in correct order
        gff_df = self.merged_df[
            [
                "contig",
                "Method",
                "Region",
                "start",
                "stop",
                "score",
                "frame",
                "phase",
                "attributes",
            ]
        ]

        # change start and stop to int
        gff_df["start"] = gff_df["start"].astype("int")
        gff_df["stop"] = gff_df["stop"].astype("int")

        ### trnas
        # check if no trnas

        # to make sure you aren't skipping trnas
        if self.skip_extra_annotations is False:
            col_list = [
                "contig",
                "Method",
                "Region",
                "start",
                "stop",
                "score",
                "frame",
                "phase",
                "attributes",
            ]
            trna_empty = is_trna_empty(self.out_dir)
            if trna_empty == False:
                trna_df = pd.read_csv(
                    os.path.join(self.out_dir, "trnascan_out.gff"),
                    delimiter="\t",
                    index_col=False,
                    names=col_list,
                )
                trna_df["contig"] = trna_df["contig"].astype(str)

                # convert the method to update with version
                trna_df["Method"] = f"tRNAscan-SE_{self.trna_version}"

                # index hack if meta mode
                if self.meta_mode == True:
                    subset_dfs = []
                    for contig in contigs:
                        subset_df = trna_df[trna_df["contig"] == contig].reset_index()
                        # keep only trnas before indexing
                        subset_df = subset_df[
                            (subset_df["Region"] == "tRNA")
                            | (subset_df["Region"] == "pseudogene")
                        ]
                        subset_df = subset_df.reset_index(drop=True)
                        subset_df["count"] = subset_df.index
                        # so not 0 indexed
                        subset_df["count"] = subset_df["count"] + 1
                        # z fill to make the locus tag 4
                        subset_df["locus_tag"] = (
                            contig
                            + "_tRNA_"
                            + subset_df["count"].astype(str).str.zfill(4)
                        )
                        subset_df = subset_df.drop(columns=["count"])
                        subset_dfs.append(subset_df)
                    trna_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
                    trna_df = trna_df.drop(columns=["index"])
                else:
                    # keep only trnas
                    trna_df = trna_df[
                        (trna_df["Region"] == "tRNA")
                        | (trna_df["Region"] == "pseudogene")
                    ]
                    trna_df = trna_df.reset_index(drop=True)
                    trna_df["count"] = trna_df.index
                    trna_df["count"] = trna_df["count"] + 1
                    trna_df["locus_tag"] = (
                        self.locustag
                        + "_tRNA_"
                        + trna_df["count"].astype(str).str.zfill(4)
                    )
                    trna_df = trna_df.drop(columns=["count"])

                trna_df.start = trna_df.start.astype(int)
                trna_df.stop = trna_df.stop.astype(int)
                trna_df[["attributes", "isotypes"]] = trna_df["attributes"].str.split(
                    ";isotype=", expand=True
                )
                trna_df[["isotypes", "anticodon"]] = trna_df["isotypes"].str.split(
                    ";anticodon=", expand=True
                )
                trna_df[["anticodon", "rest"]] = trna_df["anticodon"].str.split(
                    ";gene_biotype", expand=True
                )
                trna_df["trna_product"] = (
                    "tRNA-" + trna_df["isotypes"] + "(" + trna_df["anticodon"] + ")"
                )
                trna_df = trna_df.drop(columns=["attributes"])
                trna_df["attributes"] = (
                    "ID="
                    + trna_df["locus_tag"]
                    + ";"
                    + "transl_table="
                    + locus_df["transl_table"].astype(str)
                    + ";"
                    + "trna="
                    + trna_df["trna_product"].astype(str)
                    + ";"
                    + "isotype="
                    + trna_df["isotypes"].astype(str)
                    + ";"
                    + "anticodon="
                    + trna_df["anticodon"].astype(str)
                    + ";"
                    + "locus_tag="
                    + trna_df["locus_tag"]
                )
                trna_df = trna_df.drop(
                    columns=[
                        "isotypes",
                        "anticodon",
                        "rest",
                        "trna_product",
                        "locus_tag",
                    ]
                )

            ### crisprs
            crispr_count = get_crispr_count(self.out_dir, self.prefix)
            # add to gff if > 0
            if crispr_count > 0:
                minced_df = pd.read_csv(
                    os.path.join(self.out_dir, self.prefix + "_minced.gff"),
                    delimiter="\t",
                    index_col=False,
                    names=col_list,
                    comment="#",
                )
                minced_df["contig"] = minced_df["contig"].astype(str)
                minced_df.start = minced_df.start.astype(int)
                minced_df.stop = minced_df.stop.astype(int)
                minced_df[["attributes", "rpt_unit_seq"]] = minced_df[
                    "attributes"
                ].str.split(";rpt_unit_seq=", expand=True)
                minced_df[["attributes", "rpt_family"]] = minced_df[
                    "attributes"
                ].str.split(";rpt_family=", expand=True)
                minced_df[["attributes", "rpt_type"]] = minced_df[
                    "attributes"
                ].str.split(";rpt_type=", expand=True)
                minced_df = minced_df.drop(columns=["attributes"])
                # index hack if meta mode
                subset_dfs = []
                if self.meta_mode == True:
                    for contig in contigs:
                        subset_df = minced_df[
                            minced_df["contig"] == contig
                        ].reset_index()
                        subset_df["count"] = subset_df.index
                        # so not 0 indexed
                        subset_df["count"] = subset_df["count"] + 1
                        # z fill to make the locus tag 4
                        subset_df["count"] = subset_df["count"].astype(str).str.zfill(4)
                        subset_df["locus_tag"] = (
                            contig
                            + "_CRISPR_"
                            + subset_df["count"].astype(str).str.zfill(4)
                        )
                        subset_df = subset_df.drop(columns=["count"])
                        subset_dfs.append(subset_df)
                    minced_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
                    minced_df = minced_df.drop(columns=["index"])
                else:
                    minced_df["count"] = minced_df.index
                    minced_df["count"] = minced_df["count"] + 1
                    minced_df["locus_tag"] = (
                        self.locustag
                        + "_CRISPR_"
                        + minced_df["count"].astype(str).str.zfill(4)
                    )
                    minced_df = minced_df.drop(columns=["count"])

                minced_df["attributes"] = (
                    "ID="
                    + minced_df["locus_tag"]
                    + ";"
                    + "transl_table="
                    + locus_df["transl_table"].astype(str)
                    + ";"
                    + "rpt_type="
                    + minced_df["rpt_type"].astype(str)
                    + ";"
                    + "rpt_family="
                    + minced_df["rpt_family"].astype(str)
                    + ";"
                    + "rpt_unit_seq="
                    + minced_df["rpt_unit_seq"].astype(str)
                    + ";"
                    + "locus_tag="
                    + minced_df["locus_tag"]
                )
                minced_df = minced_df.drop(
                    columns=["rpt_unit_seq", "rpt_family", "rpt_type", "locus_tag"]
                )

            ### tmrna
            # add to gff there is a tmrna
            if self.tmrna_flag == True:
                tmrna_df = pd.read_csv(
                    os.path.join(self.out_dir, self.prefix + "_aragorn.gff"),
                    delimiter="\t",
                    index_col=False,
                    names=col_list,
                )
                tmrna_df["contig"] = tmrna_df["contig"].astype(str)
                tmrna_df.start = tmrna_df.start.astype(int)
                tmrna_df.stop = tmrna_df.stop.astype(int)

                # index hack if meta mode
                subset_dfs = []
                if self.meta_mode == True:
                    for contig in contigs:
                        subset_df = tmrna_df[tmrna_df["contig"] == contig].reset_index()
                        subset_df["count"] = subset_df.index
                        # so not 0 indexed
                        subset_df["count"] = subset_df["count"] + 1
                        # z fill to make the locus tag 4
                        subset_df["count"] = subset_df["count"].astype(str).str.zfill(4)
                        subset_df["locus_tag"] = (
                            contig
                            + "_tmRNA_"
                            + subset_df["count"].astype(str).str.zfill(4)
                        )
                        subset_df = subset_df.drop(columns=["count"])
                        subset_dfs.append(subset_df)
                    tmrna_df = pd.concat(subset_dfs, axis=0, ignore_index=True)
                    tmrna_df = tmrna_df.drop(columns=["index"])
                else:
                    tmrna_df["count"] = tmrna_df.index
                    tmrna_df["count"] = tmrna_df["count"] + 1
                    tmrna_df["locus_tag"] = (
                        self.locustag
                        + "_tmRNA_"
                        + tmrna_df["count"].astype(str).str.zfill(4)
                    )
                    tmrna_df = tmrna_df.drop(columns=["count"])

                tmrna_df["attributes"] = (
                    "ID="
                    + tmrna_df["locus_tag"]
                    + ";"
                    + "transl_table="
                    + locus_df["transl_table"].astype(str)
                    + ";"
                    + tmrna_df["attributes"].astype(str)
                    + ";locus_tag="
                    + tmrna_df["locus_tag"]
                )
                tmrna_df = tmrna_df.drop(columns=["locus_tag"])

        # write header of final gff files
        with open(os.path.join(self.out_dir, self.prefix + ".gff"), "w") as f:
            f.write("##gff-version 3\n")
            for index, row in self.length_df.iterrows():
                f.write(
                    "##sequence-region "
                    + row["contig"]
                    + " 1 "
                    + str(row["length"])
                    + "\n"
                )

        # combine dfs depending on whether the elements were detected

        # skip extra trna setc
        if self.skip_extra_annotations is True:
            df_list = [gff_df]
        else:
            if (
                trna_empty is True and self.tmrna_flag is False and crispr_count == 0
            ):  # all
                df_list = [gff_df]
            elif trna_empty is False and self.tmrna_flag is False and crispr_count == 0:
                df_list = [gff_df, trna_df]
            elif trna_empty is True and self.tmrna_flag is True and crispr_count == 0:
                df_list = [gff_df, tmrna_df]
            elif trna_empty is True and self.tmrna_flag is False and crispr_count > 0:
                df_list = [gff_df, minced_df]
            elif trna_empty is False and self.tmrna_flag is True and crispr_count == 0:
                df_list = [gff_df, trna_df, tmrna_df]
            elif trna_empty is False and self.tmrna_flag is False and crispr_count > 0:
                df_list = [gff_df, trna_df, minced_df]
            elif trna_empty is True and self.tmrna_flag is True and crispr_count > 0:
                df_list = [gff_df, tmrna_df, minced_df]
            # if trna_empty is False and self.tmrna_flag is True and crispr_count > 0:  # all detected
            else:  # all detected
                df_list = [gff_df, trna_df, tmrna_df, minced_df]

        total_gff = pd.concat(df_list, ignore_index=True)

        # ensure that the start and stops are integer
        total_gff.start = total_gff.start.astype(int)
        total_gff.stop = total_gff.stop.astype(int)

        # sorts all features by start
        total_gff = total_gff.groupby(
            ["contig"], sort=False, as_index=False, group_keys=True
        ).apply(pd.DataFrame.sort_values, "start", ascending=True)

        # write final gff to file
        with open(os.path.join(self.out_dir, self.prefix + ".gff"), "a") as f:
            total_gff.to_csv(f, sep="\t", index=False, header=False)

        # write fasta on the end
        with open(os.path.join(self.out_dir, self.prefix + ".gff"), "a") as f:
            f.write("##FASTA\n")
            fasta_sequences = SeqIO.parse(open(self.input_fasta), "fasta")
            # SeqIO.write(fasta_sequences, f, "fasta")
            for record in fasta_sequences:
                # to write the contig id not the contig description - https://github.com/gbouras13/pharokka/issues/267
                f.write(f">{record.id}\n")
                sequence = record.seq
                chunk_size = 60
                for i in range(0, len(sequence), chunk_size):
                    f.write(str(sequence[i : i + chunk_size]) + "\n")

        self.locus_df = locus_df
        self.gff_df = gff_df
        self.total_gff = total_gff
        if self.skip_extra_annotations is False:
            self.trna_empty = trna_empty
            self.crispr_count = crispr_count
        else:  # skip annotations
            self.trna_empty = True
            self.crispr_count = 0
            self.tmrna_flag = False

    def create_tbl(
        self,
    ):
        """
        Creates the pharokka.tbl file
        :param merged_df: a pandas df with the merged output
        :param length_df: a pandas df as output from get_contig_name_lengths()
        :param out_dir: output directory path
        :param prefix: output prefix
        :param gene_predictor: phanotate or prodigal
        :param tmrna_flag: output from create_gff()
        :param gff_df: dataframe for the gff creation
        :coding_table: int with the prodigal coding table used
        :return:
        """

        # get the cds

        self.total_gff = self.total_gff.reset_index(drop=True)

        if self.gene_predictor == "phanotate":
            cds_df = self.total_gff[
                self.total_gff["Method"] == f"PHANOTATE_{self.phanotate_version}"
            ]
        elif self.gene_predictor == "prodigal":
            cds_df = self.total_gff[
                self.total_gff["Method"] == f"Pyrodigal_{self.pyrodigal_version}"
            ]
        elif self.gene_predictor == "prodigal-gv":
            cds_df = self.total_gff[
                self.total_gff["Method"] == f"Pyrodigal-gv_{self.pyrodigal_gv_version}"
            ]
        elif self.gene_predictor == "genbank":
            cds_df = self.total_gff[self.total_gff["Method"] == "CUSTOM"]

        cds_df[["attributes", "locus_tag"]] = cds_df["attributes"].str.split(
            ";locus_tag=", expand=True
        )
        cds_df[["locus_tag", "rest"]] = cds_df["locus_tag"].str.split(
            ";function=", expand=True
        )

        ### trnas
        # check if no trnas
        if self.trna_empty is False:
            trna_df = self.total_gff[
                self.total_gff["Method"] == f"tRNAscan-SE_{self.trna_version}"
            ]
            # keep only trnas and pseudogenes
            trna_df.contig = trna_df.contig.astype(str)
            trna_df.start = trna_df.start.astype(int)
            trna_df.stop = trna_df.stop.astype(int)
            trna_df[["attributes", "locus_tag"]] = trna_df["attributes"].str.split(
                ";locus_tag=", expand=True
            )
            trna_df[["attributes", "isotypes"]] = trna_df["attributes"].str.split(
                ";isotype=", expand=True
            )
            trna_df[["isotypes", "anticodon"]] = trna_df["isotypes"].str.split(
                ";anticodon=", expand=True
            )
            trna_df["trna_product"] = (
                "tRNA-" + trna_df["isotypes"] + "(" + trna_df["anticodon"] + ")"
            )

        #### CRISPRs
        if self.crispr_count > 0:
            crispr_df = self.total_gff[self.total_gff["Region"] == "repeat_region"]
            crispr_df.contig = crispr_df.contig.astype(str)
            crispr_df.start = crispr_df.start.astype(int)
            crispr_df.stop = crispr_df.stop.astype(int)
            crispr_df[["attributes", "locus_tag"]] = crispr_df["attributes"].str.split(
                ";locus_tag=", expand=True
            )
            crispr_df[["attributes", "rpt_unit_seq"]] = crispr_df[
                "attributes"
            ].str.split(";rpt_unit_seq=", expand=True)

        ### TMRNAs
        if self.tmrna_flag is True:
            tmrna_df = self.total_gff[self.total_gff["Region"] == "tmRNA"]
            tmrna_df.contig = tmrna_df.contig.astype(str)
            tmrna_df.start = tmrna_df.start.astype(int)
            tmrna_df.stop = tmrna_df.stop.astype(int)
            tmrna_df[["attributes", "locus_tag"]] = tmrna_df["attributes"].str.split(
                ";locus_tag=", expand=True
            )

        with open(os.path.join(self.out_dir, self.prefix + ".tbl"), "w") as f:
            for index, row in self.length_df.iterrows():
                contig = str(row["contig"])
                f.write(">Feature " + contig + "\n")
                subset_df = self.merged_df[self.merged_df["contig"] == contig]
                for index, row in subset_df.iterrows():
                    start = str(row["start"])
                    stop = str(row["stop"])
                    if row["frame"] == "-":
                        start = str(row["stop"])
                        stop = str(row["start"])
                    f.write(start + "\t" + stop + "\t" + row["Region"] + "\n")
                    f.write(
                        ""
                        + "\t"
                        + ""
                        + "\t"
                        + ""
                        + "\t"
                        + "product"
                        + "\t"
                        + str(row["annot"])
                        + "\n"
                    )
                    f.write(
                        ""
                        + "\t"
                        + ""
                        + "\t"
                        + ""
                        + "\t"
                        + "locus_tag"
                        + "\t"
                        + str(row["locus_tag"])
                        + "\n"
                    )
                    f.write(
                        ""
                        + "\t"
                        + ""
                        + "\t"
                        + ""
                        + "\t"
                        + "transl_table"
                        + "\t"
                        + str(row["transl_table"])
                        + "\n"
                    )
                if self.trna_empty == False:
                    subset_trna_df = trna_df[trna_df["contig"] == contig]
                    for index, row in subset_trna_df.iterrows():
                        start = str(row["start"])
                        stop = str(row["stop"])
                        if row["frame"] == "-":
                            start = str(row["stop"])
                            stop = str(row["start"])
                        f.write(start + "\t" + stop + "\t" + row["Region"] + "\n")
                        f.write(
                            ""
                            + "\t"
                            + ""
                            + "\t"
                            + ""
                            + "\t"
                            + "product"
                            + "\t"
                            + str(row["trna_product"])
                            + "\n"
                        )
                        f.write(
                            ""
                            + "\t"
                            + ""
                            + "\t"
                            + ""
                            + "\t"
                            + "locus_tag"
                            + "\t"
                            + str(row["locus_tag"])
                            + "\n"
                        )
                if self.crispr_count > 0:
                    subset_crispr_df = crispr_df[crispr_df["contig"] == contig]
                    for index, row in subset_crispr_df.iterrows():
                        start = str(row["start"])
                        stop = str(row["stop"])
                        if row["frame"] == "-":
                            start = str(row["stop"])
                            stop = str(row["start"])
                        f.write(start + "\t" + stop + "\t" + row["Region"] + "\n")
                        f.write(
                            ""
                            + "\t"
                            + ""
                            + "\t"
                            + ""
                            + "\t"
                            + "locus_tag"
                            + "\t"
                            + str(row["locus_tag"])
                            + "\n"
                        )
                        f.write(
                            ""
                            + "\t"
                            + ""
                            + "\t"
                            + ""
                            + "\t"
                            + "product"
                            + "\t"
                            + str(row["rpt_unit_seq"])
                            + "\n"
                        )
                if self.tmrna_flag == True:
                    subset_tmrna_df = tmrna_df[tmrna_df["contig"] == contig]
                    for index, row in subset_tmrna_df.iterrows():
                        start = str(row["start"])
                        stop = str(row["stop"])
                        if row["frame"] == "-":
                            start = str(row["stop"])
                            stop = str(row["start"])
                        f.write(start + "\t" + stop + "\t" + "tmRNA" + "\n")
                        f.write(
                            ""
                            + "\t"
                            + ""
                            + "\t"
                            + ""
                            + "\t"
                            + "locus_tag"
                            + "\t"
                            + str(row["locus_tag"])
                            + "\n"
                        )
                        f.write(
                            ""
                            + "\t"
                            + ""
                            + "\t"
                            + ""
                            + "\t"
                            + "product"
                            + "\t"
                            + "transfer-messenger RNA, SsrA"
                            + "\n"
                        )

    def create_gff_singles(self):
        """
        Creates the single gff3 files for each inout contig in meta
        :param length_df: a pandas df as output from get_contig_name_lengths()
        :param fasta_input: input fasta file
        :param out_dir: output directory path
        :param locus_df: output from create_gff
        :param total_gff: output from create_gff
        """

        single_gff_dir = os.path.join(self.out_dir, "single_gffs")

        check_and_create_directory(single_gff_dir)

        for index, row in self.length_df.iterrows():
            contig = row["contig"]
            with open(os.path.join(single_gff_dir, contig + ".gff"), "w") as f:
                # write header of final gff files
                f.write("##gff-version 3\n")
                f.write(
                    "##sequence-region " + contig + " 1 " + str(row["length"]) + "\n"
                )

            subset_df = self.total_gff[self.total_gff["contig"] == contig]

            # write final gff to file
            with open(os.path.join(single_gff_dir, contig + ".gff"), "a") as f:
                subset_df.to_csv(f, sep="\t", index=False, header=False)

            # write fasta on the end
            with open(os.path.join(single_gff_dir, contig + ".gff"), "a") as f:
                f.write("##FASTA\n")
                fasta_sequences = SeqIO.parse(open(self.input_fasta), "fasta")
                # match id to contig header
                for dna_record in fasta_sequences:
                    if dna_record.id == contig:
                        SeqIO.write(dna_record, f, "fasta")

    def convert_singles_gff_to_gbk(self):
        """
        Converts all single gffs to gbks
        :param length_df: a pandas df as output from get_contig_name_lengths()
        :param out_dir: output directory path
        :param num_fastas: int, number of fasta in input
        """

        single_gff_dir = os.path.join(self.out_dir, "single_gffs")
        single_gbk_dir = os.path.join(self.out_dir, "single_gbks")

        check_and_create_directory(single_gbk_dir)

        # directory of all split fastas
        split_fasta_dir = os.path.join(self.out_dir, "input_split_tmp")

        for index, row in self.length_df.iterrows():
            # get the input fasta for each contig
            fasta_file = os.path.join(
                split_fasta_dir, "input_subprocess" + str(index + 1) + ".fasta"
            )
            contig = row["contig"]
            convert_gff_to_gbk(
                fasta_file, single_gff_dir, single_gbk_dir, contig, self.prot_seq_df
            )

    def split_fasta_singles(self):
        """Splits the input fasta into separate single fasta files for output based on contig names

        :param input_fasta: input multifasta file
        :param out_dir: output directory
        """

        single_fastas = os.path.join(self.out_dir, "single_fastas")
        check_and_create_directory(single_fastas)
        fasta_sequences = SeqIO.parse(open(self.input_fasta), "fasta")

        for dna_record in fasta_sequences:
            contig = dna_record.id
            with open(os.path.join(single_fastas, contig + ".fasta"), "w") as f:
                SeqIO.write(dna_record, f, "fasta")

    def split_faas_singles(self):
        """Splits the .faa fasta into separate single fasta files for output based on contig names

        :param input_fasta: input multifasta file
        :param out_dir: output directory
        """

        single_faas = os.path.join(self.out_dir, "single_faas")
        check_and_create_directory(single_faas)
        faa_file = os.path.join(self.out_dir, f"{self.gene_predictor}.faa")
        faa_sequences = SeqIO.parse(open(faa_file), "fasta")

        for record in faa_sequences:
            # remove the last 9 chars e.g. _CDS_0001
            protein_id = record.id[:-9]
            # needs to be append
            with open(os.path.join(single_faas, f"{protein_id}.faa"), "a") as f:
                SeqIO.write(record, f, "fasta")

    def write_tophits_vfdb_card(self):
        """
        Outputs top_hits_vfdb.tsv and top_hits_card.tsv
        :param merged_df: a pandas df as output from process_results()
        :param locus_df: a pandas df as output from create_gff()
        :param out_dir: output directory path
        :param prefix: output prefix
        :return:
        """
        ######################################
        ##### update vfdb with locus tag #####
        ###### merge locus into the vfdb ####
        ######################################
        # needed for vfdb card matching later

        genes_for_vfdb_card = self.merged_df["gene"]
        self.locus_df["gene"] = genes_for_vfdb_card
        self.vfdb_results = self.vfdb_results.merge(
            self.locus_df, how="left", on="gene"
        )

        # get a list of columns
        cols = list(self.vfdb_results)

        # move the column to head of list using index, pop and insert
        cols.insert(0, cols.pop(cols.index("locus_tag")))
        self.vfdb_results = self.vfdb_results.loc[:, cols]

        # keep only desired columns  and save
        self.vfdb_results = self.vfdb_results[
            [
                "contig",
                "locus_tag",
                "vfdb_hit_x",
                "vfdb_alnScore_x",
                "vfdb_seqIdentity_x",
                "start",
                "stop",
                "frame",
            ]
        ]

        self.vfdb_results.columns = [
            "contig",
            "gene",
            "vfdb_hit",
            "vfdb_alnScore",
            "vfdb_seqIdentity",
            "start",
            "stop",
            "frame",
        ]

        if self.mmseqs_flag is True:
            self.vfdb_results = self.vfdb_results.sort_values(by=["start"])
            self.vfdb_results.to_csv(
                os.path.join(self.out_dir, "top_hits_vfdb.tsv"), sep="\t", index=False
            )

        ######################################
        ##### update card with locus tag #####
        ###### merge locus into the card df ####
        #####################################
        self.card_results = self.card_results.merge(
            self.locus_df, how="left", on="gene"
        )
        # get a list of columns
        cols = list(self.card_results)
        # move the column to head of list using index, pop and insert
        cols.insert(0, cols.pop(cols.index("locus_tag")))
        self.card_results = self.card_results.loc[:, cols]

        # keep only desired columns   sand save
        self.card_results = self.card_results[
            [
                "contig",
                "locus_tag",
                "CARD_hit_x",
                "CARD_alnScore_x",
                "CARD_seqIdentity_x",
                "start",
                "stop",
                "frame",
            ]
        ]
        self.card_results.columns = [
            "contig",
            "gene",
            "card_hit",
            "card_alnScore",
            "card_seqIdentity",
            "start",
            "stop",
            "frame",
        ]
        if self.mmseqs_flag is True:
            self.card_results = self.card_results.sort_values(by=["start"])
            self.card_results.to_csv(
                os.path.join(self.out_dir, "top_hits_card.tsv"), sep="\t", index=False
            )

    def create_txt(self):
        """
        Creates the _cds_functions.tsv and _length_gc_cds_density.tsv outputs
        :param merged_df: a pandas df as output from process_results()
        :param length_df: a pandas df as output from get_contig_name_lengths()
        :param out_dir: output directory path
        :param prefix: output prefix
        :return:
        """

        # get contigs - convert to string to make the match work for integer contigs
        contigs = self.length_df["contig"].astype("string")
        # convert to string to make the match work for integer contigs
        self.merged_df["contig"] = self.merged_df["contig"].astype("string")

        # list of all dataframes with functions (for later)
        combo_list = []

        #### trna scan
        # read in trnascan
        if self.skip_extra_annotations is False:
            col_list = [
                "contig",
                "Method",
                "Region",
                "start",
                "stop",
                "score",
                "frame",
                "phase",
                "attributes",
            ]
            trna_df = pd.read_csv(
                os.path.join(self.out_dir, "trnascan_out.gff"),
                delimiter="\t",
                index_col=False,
                names=col_list,
            )
            # keep only trnas and pseudogenes
            trna_df = trna_df[
                (trna_df["Region"] == "tRNA") | (trna_df["Region"] == "pseudogene")
            ]

            #### crispr
            crispr_df = pd.read_csv(
                os.path.join(self.out_dir, self.prefix + "_minced.gff"),
                delimiter="\t",
                index_col=False,
                names=col_list,
                comment="#",
            )

            #### tmrna
            tmrna_df = pd.read_csv(
                os.path.join(self.out_dir, self.prefix + "_aragorn.gff"),
                delimiter="\t",
                index_col=False,
                names=col_list,
            )

        # write descriptions for each contig
        for contig in contigs:
            # get cds's in the contig
            cds_mmseqs_merge_cont_df = self.merged_df[
                self.merged_df["contig"] == contig
            ]
            # counts of the cds trnas
            cds_count = len(
                cds_mmseqs_merge_cont_df[cds_mmseqs_merge_cont_df["Region"] == "CDS"]
            )
            if self.skip_extra_annotations is False:
                trna_count = len(trna_df[trna_df["contig"] == contig])
                tmrna_count = len(tmrna_df[tmrna_df["contig"] == contig])
                crispr_count = len(crispr_df[crispr_df["contig"] == contig])
            if len(self.vfdb_results["contig"]) != 0:
                vfdb_count = len(
                    self.vfdb_results[self.vfdb_results["contig"] == contig]
                )
                vfdb_count = self.vfdb_results["contig"].str.contains(contig).sum()
            else:
                vfdb_count = 0
            if len(self.card_results["contig"]) != 0:
                CARD_count = len(
                    self.card_results[self.card_results["contig"] == contig]
                )
                CARD_count = self.card_results["contig"].str.contains(contig).sum()
            else:
                CARD_count = 0
            # get the total length of the contig
            contig_length = self.length_df[self.length_df["contig"] == contig]["length"]
            if cds_count > 0:
                # gets the total cds coding length
                cds_lengths = abs(
                    cds_mmseqs_merge_cont_df["start"] - cds_mmseqs_merge_cont_df["stop"]
                ).sum()
                # get function
                cds_mmseqs_merge_cont_df[["attributes2"]] = cds_mmseqs_merge_cont_df[
                    ["attributes"]
                ]
                cds_mmseqs_merge_cont_df[
                    ["attributes2", "function"]
                ] = cds_mmseqs_merge_cont_df["attributes2"].str.split(
                    ";function=", expand=True
                )
                cds_mmseqs_merge_cont_df = cds_mmseqs_merge_cont_df.drop(
                    columns=["attributes2"]
                )
                cds_mmseqs_merge_cont_df[
                    ["function", "product"]
                ] = cds_mmseqs_merge_cont_df["function"].str.split(
                    ";product=", expand=True
                )
                cds_mmseqs_merge_cont_df = cds_mmseqs_merge_cont_df.drop(
                    columns=["product"]
                )

                # get counts of functions and cds
                # all 10 PHROGs categories
                # modified for v1.0.0 in case some have 0s (to make consistent for downstream)
                connector_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"] == "connector"
                    ]
                )
                metabolism_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"]
                        == "DNA, RNA and nucleotide metabolism"
                    ]
                )
                head_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"] == "head and packaging"
                    ]
                )
                integration_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"]
                        == "integration and excision"
                    ]
                )
                lysis_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"] == "lysis"
                    ]
                )
                moron_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"]
                        == "moron, auxiliary metabolic gene and host takeover"
                    ]
                )
                other_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"] == "other"
                    ]
                )
                tail_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"] == "tail"
                    ]
                )
                transcription_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"]
                        == "transcription regulation"
                    ]
                )
                unknown_count = len(
                    cds_mmseqs_merge_cont_df[
                        cds_mmseqs_merge_cont_df["function"] == "unknown function"
                    ]
                )
                # create count list  for the dataframe
                count_list = [
                    cds_count,
                    connector_count,
                    metabolism_count,
                    head_count,
                    integration_count,
                    lysis_count,
                    moron_count,
                    other_count,
                    tail_count,
                    transcription_count,
                    unknown_count,
                ]
            else:
                cds_lengths = 0

                count_list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

            description_list = [
                "CDS",
                "connector",
                "DNA, RNA and nucleotide metabolism",
                "head and packaging",
                "integration and excision",
                "lysis",
                "moron, auxiliary metabolic gene and host takeover",
                "other",
                "tail",
                "transcription regulation",
                "unknown function",
            ]
            contig_list = [
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
                contig,
            ]
            # cds df
            cds_df = pd.DataFrame(
                {
                    "Description": description_list,
                    "Count": count_list,
                    "contig": contig_list,
                }
            )

            # only if not skipped
            if self.skip_extra_annotations is False:
                # add other features
                trna_row = pd.DataFrame(
                    {
                        "Description": ["tRNAs"],
                        "Count": [trna_count],
                        "contig": [contig],
                    }
                )
                crispr_row = pd.DataFrame(
                    {
                        "Description": ["CRISPRs"],
                        "Count": [crispr_count],
                        "contig": [contig],
                    }
                )
                tmrna_row = pd.DataFrame(
                    {
                        "Description": ["tmRNAs"],
                        "Count": [tmrna_count],
                        "contig": [contig],
                    }
                )

            vfdb_row = pd.DataFrame(
                {
                    "Description": ["VFDB_Virulence_Factors"],
                    "Count": [vfdb_count],
                    "contig": [contig],
                }
            )
            CARD_row = pd.DataFrame(
                {
                    "Description": ["CARD_AMR_Genes"],
                    "Count": [CARD_count],
                    "contig": [contig],
                }
            )
            # calculate the cds coding density and add to length_df
            cds_coding_density = cds_lengths * 100 / contig_length
            cds_coding_density = round(cds_coding_density, 2)
            self.length_df.loc[
                self.length_df["contig"] == contig, "cds_coding_density"
            ] = cds_coding_density
            # eappend it all to combo_list
            combo_list.append(cds_df)
            if self.skip_extra_annotations is False:
                combo_list.append(trna_row)
                combo_list.append(crispr_row)
                combo_list.append(tmrna_row)
            combo_list.append(vfdb_row)
            combo_list.append(CARD_row)

        # combine all contigs into one final df
        description_total_df = pd.concat(combo_list)
        # save the output
        description_total_df.to_csv(
            os.path.join(self.out_dir, self.prefix + "_cds_functions.tsv"),
            sep="\t",
            index=False,
        )
        # save the length_gc.tsv also
        self.length_df.to_csv(
            os.path.join(self.out_dir, self.prefix + "_length_gc_cds_density.tsv"),
            sep="\t",
            index=False,
        )

    def update_fasta_headers(self):
        """
        Updates the fasta output headers to have a consistent locus tag & gene description for downstrea use
        :param locus_df a pandas df as output from create_gff()
        :param out_dir: output directory path
        :gene_predictor: string 'phanotate' or 'prodigal' with the gene predictor used
        """

        # define outputs
        fasta_input_nts_tmp = self.gene_predictor + "_out_tmp.fasta"
        fasta_input_aas_tmp = self.gene_predictor + "_aas_tmp.fasta"
        fasta_output_nts_gd = self.gene_predictor + ".ffn"
        fasta_output_aas_gd = self.gene_predictor + ".faa"

        # nucleotides

        with open(os.path.join(self.out_dir, fasta_output_nts_gd), "w") as nt_fa:
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_nts_tmp), "fasta"
            ):
                dna_record.id = str(self.locus_df["locus_tag"].iloc[i])
                dna_record.description = str(self.locus_df["annot"].iloc[i])
                SeqIO.write(dna_record, nt_fa, "fasta")
                i += 1

        # amino acids

        with open(os.path.join(self.out_dir, fasta_output_aas_gd), "w") as aa_fa:
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_aas_tmp), "fasta"
            ):
                dna_record.id = str(self.locus_df["locus_tag"].iloc[i])
                dna_record.description = str(self.locus_df["annot"].iloc[i])
                SeqIO.write(dna_record, aa_fa, "fasta")
                i += 1

    def update_final_output(self):
        """
        Updates the fasta output headers to have a consistent locus tag & gene description for downstrea use
        :param merged_df: a pandas df as output from process_results()
        :param locus_df: a pandas df as output from create_gff()
        :param out_dir: output directory path
        :param prefix: output prefix
        :return:
        """
        # return back the merged_df but with the locus tag instead of gene
        # rename gene with locus_tag
        locus_tag_series = self.locus_df["locus_tag"]
        self.merged_df["gene"] = locus_tag_series

        #########
        # rearrange start and stop for neg strang
        #########

        st_cols = ["start", "stop"]
        # indices where start is greater than stop
        ixs = self.merged_df["frame"] == "-"
        # Where ixs is True, values are swapped
        self.merged_df.loc[ixs, st_cols] = (
            self.merged_df.loc[ixs, st_cols].reindex(columns=st_cols[::-1]).values
        )

        # get a list of columns
        cols = list(self.merged_df)
        # move the column to head of list using index, pop and insert
        cols.insert(0, cols.pop(cols.index("gene")))
        self.merged_df = self.merged_df.loc[:, cols]
        # drop cols
        self.merged_df = self.merged_df.drop(
            columns=["phase", "attributes", "count", "locus_tag"]
        )

        # write output
        final_output_path = os.path.join(
            self.out_dir, self.prefix + "_cds_final_merged_output.tsv"
        )
        self.merged_df.to_csv(final_output_path, sep="\t", index=False)

    def extract_terl(self):
        """
        Extract large terminase subunit
        :param locus_df a pandas df as output from create_gff()
        :param out_dir: output directory path
        :gene_predictor: string 'phanotate' or 'prodigal' with the gene predictor used
        """

        # phanotate
        fasta_input_nts_tmp = self.gene_predictor + "_out_tmp.fasta"
        fasta_input_aas_tmp = self.gene_predictor + "_aas_tmp.fasta"

        # nucleotide

        with open(os.path.join(self.out_dir, "terL.ffn"), "w") as aa_fa:
            # i loops over the rows of the dataframe
            # j counts the number of terLs
            i = 0
            j = 0
            for dna_record in SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_nts_tmp), "fasta"
            ):
                dna_record.id = str(self.locus_df["locus_tag"].iloc[i])
                dna_record.description = str(self.locus_df["annot"].iloc[i])
                if self.locus_df["annot"].iloc[i] == "terminase large subunit":
                    SeqIO.write(dna_record, aa_fa, "fasta")
                    # report terL found the first time
                    if j < 1:
                        logger.info("Terminase large subunit found.")
                        # report multiple found the second time
                    if j == 1:
                        logger.info(
                            "More than one CDS annotated as terminase large subunit found. \nSaving all."
                        )
                    j += 1
                i += 1

        # amino acid no need to print

        with open(os.path.join(self.out_dir, "terL.faa"), "w") as aa_fa:
            # i loops over the rows of the dataframe
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_aas_tmp), "fasta"
            ):
                dna_record.id = str(self.locus_df["locus_tag"].iloc[i])
                dna_record.description = str(self.locus_df["annot"].iloc[i])
                if self.locus_df["annot"].iloc[i] == "terminase large subunit":
                    SeqIO.write(dna_record, aa_fa, "fasta")
                i += 1

    def inphared_top_hits(self):
        """
        Process mash output to get inphared top hits
        :param out_dir: output directory
        :param length_df: a pandas df as output from get_contig_name_lengths()
        """

        mash_tsv = os.path.join(self.out_dir, "mash_out.tsv")
        col_list = [
            "contig",
            "Accession",
            "mash_distance",
            "mash_pval",
            "mash_matching_hashes",
        ]

        # get contigs - convert to string to make the match work for integer contigs
        contigs = self.length_df["contig"].astype(str)

        # instantiate tophits list
        tophits_mash_df = []

        # read in the mash output
        mash_df = pd.read_csv(mash_tsv, delimiter="\t", index_col=False, names=col_list)
        # as string to match the contigs
        mash_df["contig"] = mash_df["contig"].astype(str)

        # instantiate tophits list
        tophits = []

        for contig in contigs:
            hit_df = (
                mash_df.loc[mash_df["contig"] == contig]
                .sort_values("mash_distance")
                .reset_index(drop=True)
            )
            hits = len(hit_df["mash_distance"])
            # add only if there is a hit
            if hits > 0:
                # top hit
                top_df = (
                    mash_df.loc[mash_df["contig"] == contig]
                    .sort_values("mash_distance")
                    .reset_index(drop=True)
                    .loc[0]
                )
                tophits.append(
                    [
                        top_df.contig,
                        top_df.Accession,
                        top_df.mash_distance,
                        top_df.mash_pval,
                        top_df.mash_matching_hashes,
                    ]
                )
            else:
                tophits.append(
                    [
                        contig,
                        "no_inphared_mash_hit",
                        "no_inphared_mash_hit",
                        "no_inphared_mash_hit",
                        "no_inphared_mash_hit",
                    ]
                )
            # create tophits df
        tophits_mash_df = pd.DataFrame(
            tophits,
            columns=[
                "contig",
                "Accession",
                "mash_distance",
                "mash_pval",
                "mash_matching_hashes",
            ],
        )

        # read in the plasdb tsv
        inphared_tsv_file = os.path.join(self.db_dir, "1Aug2023_data.tsv")

        cols = [
            "Accession",
            "Description",
            "Classification",
            "Genome_Length_(bp)",
            "Jumbophage",
            "molGC_(%)",
            "Molecule",
            "Modification_Date",
            "Number_CDS",
            "Positive_Strand_(%)",
            "Negative_Strand_(%)",
            "Coding_Capacity_(%)",
            "Low_Coding_Capacity_Warning",
            "tRNAs",
            "Host",
            "Lowest_Taxa",
            "Genus",
            "Sub-family",
            "Family",
            "Order",
            "Class",
            "Phylum",
            "Kingdom",
            "Realm",
            "Baltimore_Group",
            "Genbank_Division",
            "Isolation_Host_(beware_inconsistent_and_nonsense_values)",
        ]
        inphared_tsv_file = pd.read_csv(
            inphared_tsv_file,
            delimiter="\t",
            index_col=False,
            names=cols,
            skiprows=1,
            low_memory=False,
        )

        combined_df = tophits_mash_df.merge(
            inphared_tsv_file, on="Accession", how="left"
        )

        combined_df.to_csv(
            os.path.join(self.out_dir, self.prefix + "_top_hits_mash_inphared.tsv"),
            sep="\t",
            index=False,
        )


#########################################
####### non class functions #############
########################################


def create_mmseqs_tophits(out_dir):
    """
    creates tophits_df dataframe from mmseqs2 phrog results
    """

    ##mmseqs
    mmseqs_file = os.path.join(out_dir, "mmseqs_results.tsv")
    logger.info("Processing MMseqs2 outputs.")
    logger.info("Processing PHROGs output.")
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

    # optimise the tophits generation
    # Group by 'gene' and find the top hit for each group
    tophits_df = (
        mmseqs_df.groupby("gene", group_keys=True)
        .apply(lambda group: group.nsmallest(1, "mmseqs_eVal"))
        .reset_index(drop=True)
    )

    # create tophits df
    tophits_df = tophits_df[
        [
            "mmseqs_phrog",
            "gene",
            "mmseqs_alnScore",
            "mmseqs_seqIdentity",
            "mmseqs_eVal",
        ]
    ]

    tophits_df.to_csv(
        os.path.join(out_dir, "top_hits_mmseqs.tsv"), sep="\t", index=False
    )
    return tophits_df


def remove_post_processing_files(out_dir, gene_predictor, meta):
    """
    Cleans temporary files up
    :param out_dir: output directory path
    :param gene_predictor: phanotate or prodigal
    :return:
    """
    remove_directory(os.path.join(out_dir, "target_dir"))
    remove_directory(os.path.join(out_dir, "tmp_dir"))
    remove_directory(os.path.join(out_dir, "mmseqs"))
    remove_directory(os.path.join(out_dir, "vfdb_tmp_dir"))
    remove_directory(os.path.join(out_dir, "VFDB_dir"))
    remove_directory(os.path.join(out_dir, "VFDB"))
    remove_directory(os.path.join(out_dir, "vfdb"))
    remove_file(os.path.join(out_dir, "vfdb_results.tsv"))
    remove_directory(os.path.join(out_dir, "CARD_tmp_dir"))
    remove_directory(os.path.join(out_dir, "CARD"))
    remove_directory(os.path.join(out_dir, "CARD_dir"))
    remove_file(os.path.join(out_dir, "CARD_results.tsv"))
    remove_file(os.path.join(out_dir, "cleaned_" + gene_predictor + ".tsv"))
    remove_file(os.path.join(out_dir, "input_fasta_delim.fasta"))
    remove_file(os.path.join(out_dir, "mmseqs_results.tsv"))
    remove_file(os.path.join(out_dir, "top_hits_mmseqs.tsv"))

    # leave in tophits
    remove_file(os.path.join(out_dir, gene_predictor + "_aas_tmp.fasta"))
    remove_file(os.path.join(out_dir, gene_predictor + "_out_tmp.fasta"))
    remove_file(os.path.join(out_dir, "pharokka_tmp.gff"))
    remove_file(os.path.join(out_dir, "mash_out.tsv"))
    remove_file(os.path.join(out_dir, "input_mash_sketch.msh"))

    if gene_predictor == "phanotate":
        remove_file(os.path.join(out_dir, "phanotate_out.txt"))
    if gene_predictor == "prodigal":
        remove_file(os.path.join(out_dir, "prodigal_out.gff"))
        remove_file(os.path.join(out_dir, "prodigal_out_aas_tmp.fasta"))
    elif gene_predictor == "prodigal-gv":
        remove_file(os.path.join(out_dir, "prodigal-gv_out.gff"))
        remove_file(os.path.join(out_dir, "prodigal-gv_out_aas_tmp.fasta"))
    # delete the tmp meta files
    if meta == True:
        remove_directory(os.path.join(out_dir, "input_split_tmp/"))

    # if genbank input
    remove_file(os.path.join(out_dir, "genbank.fasta"))


def get_crispr_count(out_dir, prefix):
    """
    Gets number of crisprs
    :param out_dir: output directory path
    :param prefix: prefix
    :return: crispr_count integer
    """
    crispr_file = os.path.join(out_dir, prefix + "_minced.gff")
    with open(crispr_file) as file:
        lines = file.readlines()
    crispr_count = 0
    for line in lines:
        if line[0] != "#":
            crispr_count += 1
    return crispr_count


# check if the trna file has more than 1 line (not empty)
def is_trna_empty(out_dir):
    """
    Determines if trna output file is empty
    :param out_dir: output directory path
    :return: trna_empty Boolean
    """
    trna_empty = False
    if os.stat(os.path.join(out_dir, "trnascan_out.gff")).st_size == 0:
        trna_empty = True
    return trna_empty


#### process pyhmmer hits
def process_pyhmmer_results(merged_df, pyhmmer_results_dict):
    """
    Processes pyhmmer
    :param merged_df: merged_df in process_results
    :param pyhmmer_results_dict: dictionary with pyhmmer results
    :return: merged_df merged_df updated with pyhmmer results
    """

    # split to get protein name to match with pyhmmer
    merged_df["temp_prot"] = merged_df["gene"].str.split(" ", n=1).str[0]

    merged_df["pyhmmer_phrog"] = "No_PHROGs_HMM"
    merged_df["pyhmmer_bitscore"] = "No_PHROGs_HMM"
    merged_df["pyhmmer_evalue"] = "No_PHROGs_HMM"

    # iterate over the df
    for index, row in merged_df.iterrows():
        if (
            row["temp_prot"] in pyhmmer_results_dict
        ):  # check if the protein is in the dictionary
            merged_df.at[index, "pyhmmer_phrog"] = pyhmmer_results_dict[
                row["temp_prot"]
            ].phrog
            merged_df.at[index, "pyhmmer_bitscore"] = round(
                pyhmmer_results_dict[row["temp_prot"]].bitscore, 6
            )
            merged_df.at[index, "pyhmmer_evalue"] = pyhmmer_results_dict[
                row["temp_prot"]
            ].evalue

    # drop temp prot

    merged_df = merged_df.drop(columns=["temp_prot"])

    return merged_df


#### process pyhmmer hits
def process_custom_pyhmmer_results(merged_df, custom_pyhmmer_results_dict):
    """
    Processes pyhmmer results for custom db
    :param merged_df: merged_df in process_results
    :param pyhmmer_results_dict: dictionary with pyhmmer results
    :return: merged_df merged_df updated with pyhmmer results
    """

    # split to get protein name to match with pyhmmer
    merged_df["temp_prot"] = merged_df["gene"].str.split(" ", n=1).str[0]
    merged_df["custom_hmm_id"] = "No_custom_HMM"
    merged_df["custom_hmm_bitscore"] = "No_custom_HMM"
    merged_df["custom_hmm_evalue"] = "No_custom_HMM"

    # iterate over the df
    for index, row in merged_df.iterrows():
        if (
            row["temp_prot"] in custom_pyhmmer_results_dict
        ):  # check if the protein is in the dictionary
            merged_df.at[index, "custom_hmm_id"] = custom_pyhmmer_results_dict[
                row["temp_prot"]
            ].custom_hmm_id
            merged_df.at[index, "custom_hmm_bitscore"] = round(
                custom_pyhmmer_results_dict[row["temp_prot"]].bitscore, 6
            )
            merged_df.at[index, "custom_hmm_evalue"] = custom_pyhmmer_results_dict[
                row["temp_prot"]
            ].evalue

    # drop temp prot
    merged_df = merged_df.drop(columns=["temp_prot"])
    return merged_df


#### process vfdb files
def process_vfdb_results(out_dir, merged_df, proteins_flag=False):
    """
    Processes VFDB results
    :param out_dir: output directory path
    :param merged_df: merged_df in process_results
    :proteins_flag bool, True if pharokka_proteins is run because we need to strip off everything before the first space
    :return: merged_df merged_df updated with VFDB results
    """
    ##vfdb
    vfdb_file = os.path.join(out_dir, "vfdb_results.tsv")
    logger.info("Processing VFDB output.")
    col_list = [
        "vfdb_hit",
        "gene",
        "vfdb_alnScore",
        "vfdb_seqIdentity",
        "vfdb_eVal",
        "qStart",
        "qEnd",
        "qLen",
        "tStart",
        "tEnd",
        "tLen",
    ]

    # touch the file in case it doesn't exist (for --fast mode)
    touch_file(vfdb_file)

    vfdb_df = pd.read_csv(vfdb_file, delimiter="\t", index_col=False, names=col_list)

    # optimise the tophits generation
    # Group by 'gene' and find the top hit for each group
    tophits_df = (
        vfdb_df.groupby("gene", group_keys=True)
        .apply(lambda group: group.nsmallest(1, "vfdb_eVal"))
        .reset_index(drop=True)
    )

    # create tophits df
    tophits_df = tophits_df[
        [
            "vfdb_hit",
            "gene",
            "vfdb_alnScore",
            "vfdb_seqIdentity",
            "vfdb_eVal",
        ]
    ]

    # left join vfdb to merged_df
    tophits_df["gene"] = tophits_df["gene"].astype(str)
    # error #300 - bad merge on gene
    if proteins_flag is True:
        tophits_df["gene"] = tophits_df["gene"].str.split(" ").str.get(0)

    # merge top hits into the merged df
    merged_df = merged_df.merge(tophits_df, on="gene", how="left")
    merged_df["vfdb_hit"] = merged_df["vfdb_hit"].replace(nan, "None", regex=True)
    merged_df["vfdb_alnScore"] = merged_df["vfdb_alnScore"].replace(
        nan, "None", regex=True
    )
    merged_df["vfdb_seqIdentity"] = merged_df["vfdb_seqIdentity"].replace(
        nan, "None", regex=True
    )
    merged_df["vfdb_eVal"] = merged_df["vfdb_eVal"].replace(nan, "None", regex=True)

    # if there is a hit extract information about it
    if len(tophits_df["vfdb_hit"]) > 0:
        number_vfs = len(tophits_df["vfdb_hit"])
        logger.info(str(number_vfs) + " VFDB virulence factors identified.")
        merged_df[["genbank", "desc_tmp", "vfdb_species"]] = merged_df[
            "vfdb_hit"
        ].str.split("[", expand=True)
        merged_df["vfdb_species"] = merged_df["vfdb_species"].str.replace("]", "")
        merged_df["vfdb_species"] = merged_df["vfdb_species"].str.strip()
        # genbank has the info
        merged_df["vfdb_short_name"] = merged_df["genbank"].str.split(")", n=1).str[1]
        merged_df["vfdb_description"] = merged_df["genbank"].str.split(")", n=2).str[2]
        merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].str.replace("(", "")
        merged_df["vfdb_short_name"] = (
            merged_df["vfdb_short_name"].str.split(")", n=1).str[0]
        )
        merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].str.strip()
        merged_df["vfdb_description"] = merged_df["vfdb_description"].str.strip()
        # remove and add None
        merged_df = merged_df.drop(columns=["genbank", "desc_tmp"])
        merged_df["vfdb_short_name"] = merged_df["vfdb_short_name"].replace(
            nan, "None", regex=True
        )
        merged_df["vfdb_description"] = merged_df["vfdb_description"].replace(
            nan, "None", regex=True
        )
        merged_df["vfdb_species"] = merged_df["vfdb_species"].replace(
            nan, "None", regex=True
        )
    else:
        logger.info("0 VFDB virulence factors identified.")
        merged_df["vfdb_short_name"] = "None"
        merged_df["vfdb_description"] = "None"
        merged_df["vfdb_species"] = "None"
    return (merged_df, tophits_df)


#### process CARD files
def process_card_results(out_dir, merged_df, db_dir, proteins_flag=False):
    """
    Processes card results
    :param out_dir: output directory path
    :param merged_df: merged_df in process_results
    :proteins_flag bool, True if pharokka_proteins is run because we need to strip off everything before the first space
    :return: merged_df merged_df updated with card results
    """
    ##card
    card_file = os.path.join(out_dir, "CARD_results.tsv")
    logger.info("Processing CARD output.")
    col_list = [
        "CARD_hit",
        "gene",
        "CARD_alnScore",
        "CARD_seqIdentity",
        "CARD_eVal",
        "qStart",
        "qEnd",
        "qLen",
        "tStart",
        "tEnd",
        "tLen",
    ]
    touch_file(card_file)
    card_df = pd.read_csv(card_file, delimiter="\t", index_col=False, names=col_list)

    #
    tophits_df = (
        card_df.groupby("gene", group_keys=True)
        .apply(lambda group: group.nsmallest(1, "CARD_eVal"))
        .reset_index(drop=True)
    )

    # create tophits df
    tophits_df = tophits_df[
        [
            "CARD_hit",
            "gene",
            "CARD_alnScore",
            "CARD_seqIdentity",
            "CARD_eVal",
        ]
    ]

    # left join tophits_df to merged_df
    tophits_df["gene"] = tophits_df["gene"].astype(str)

    # error #300 - bad merge on gene in proteins, as it takes the entire FASTA headers (which may include spaces)
    # therefore take only everything before the space
    if proteins_flag is True:
        tophits_df["gene"] = tophits_df["gene"].str.split(" ").str.get(0)

    # merge top hits into the merged df
    merged_df = merged_df.merge(tophits_df, on="gene", how="left")
    merged_df["CARD_hit"] = merged_df["CARD_hit"].replace(nan, "None", regex=True)
    merged_df["CARD_alnScore"] = merged_df["CARD_alnScore"].replace(
        nan, "None", regex=True
    )
    merged_df["CARD_seqIdentity"] = merged_df["CARD_seqIdentity"].replace(
        nan, "None", regex=True
    )
    merged_df["CARD_eVal"] = merged_df["CARD_eVal"].replace(nan, "None", regex=True)
    # if there is a hit extract info
    if len(tophits_df["CARD_hit"]) > 0:
        number_cards = len(tophits_df["CARD_hit"])
        logger.info(str(number_cards) + " CARD AMR genes identified.")
        merged_df[["genbank", "CARD_species"]] = merged_df["CARD_hit"].str.split(
            "[", expand=True
        )
        merged_df["CARD_species"] = merged_df["CARD_species"].str.replace("]", "")
        merged_df["CARD_species"] = merged_df["CARD_species"].str.strip()
        merged_df[["gb", "genbank", "ARO_Accession", "CARD_short_name"]] = merged_df[
            "genbank"
        ].str.split("|", expand=True)
        merged_df["CARD_short_name"] = merged_df["CARD_short_name"].str.strip()
        # read in aro_index
        CARD_index_file = os.path.join(db_dir, "aro_index.tsv")
        col_list = [
            "ARO_Accession",
            "CVTERM_ID",
            "Model_Sequence_ID",
            "Model_ID",
            "Model_Name",
            "ARO_Name",
            "Protein_Accession",
            "DNA_Accession",
            "AMR_Gene_Family",
            "Drug_Class",
            "Resistance_Mechanism",
            "CARD_Short_Name",
        ]
        card_index_df = pd.read_csv(
            CARD_index_file, delimiter="\t", index_col=False, names=col_list, skiprows=1
        )
        card_index_df = card_index_df.drop(
            columns=[
                "CVTERM_ID",
                "Model_Sequence_ID",
                "Model_ID",
                "Model_Name",
                "ARO_Name",
                "CARD_Short_Name",
            ]
        )
        merged_df = merged_df.merge(card_index_df, on="ARO_Accession", how="left")
        merged_df = merged_df.drop(columns=["gb", "genbank"])
        merged_df["CARD_species"] = merged_df["CARD_species"].replace(
            nan, "None", regex=True
        )
        merged_df["ARO_Accession"] = merged_df["ARO_Accession"].replace(
            nan, "None", regex=True
        )
        merged_df["CARD_short_name"] = merged_df["CARD_short_name"].replace(
            nan, "None", regex=True
        )
        merged_df["Protein_Accession"] = merged_df["Protein_Accession"].replace(
            nan, "None", regex=True
        )
        merged_df["DNA_Accession"] = merged_df["DNA_Accession"].replace(
            nan, "None", regex=True
        )
        merged_df["AMR_Gene_Family"] = merged_df["AMR_Gene_Family"].replace(
            nan, "None", regex=True
        )
        merged_df["Drug_Class"] = merged_df["Drug_Class"].replace(
            nan, "None", regex=True
        )
        merged_df["Resistance_Mechanism"] = merged_df["Resistance_Mechanism"].replace(
            nan, "None", regex=True
        )
    # if no hits just Nones
    else:
        logger.info("0 CARD AMR genes identified.")
        merged_df["CARD_species"] = "None"
        merged_df["ARO_Accession"] = "None"
        merged_df["CARD_short_name"] = "None"
        merged_df["Protein_Accession"] = "None"
        merged_df["DNA_Accession"] = "None"
        merged_df["AMR_Gene_Family"] = "None"
        merged_df["Drug_Class"] = "None"
        merged_df["Resistance_Mechanism"] = "None"

    return (merged_df, tophits_df)


# check if a file is empt
def is_file_empty(file):
    """
    Determines if file is empty
    :param file: file path
    :return: empty Boolean
    """
    empty = False
    if os.stat(file).st_size == 0:
        empty = True
    return empty


def check_and_create_directory(directory):
    """Checks if directory exists, creates it if not
    :param directory: director
    """
    if os.path.isdir(directory) == False:
        os.mkdir(directory)
