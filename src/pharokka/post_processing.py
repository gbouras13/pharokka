import collections
import os
import random
import string
from pathlib import Path

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from loguru import logger

from .processes import convert_gff_to_gbk
from .util import (
    parse_attributes,
    remove_directory,
    remove_file,
    rename_file,
    touch_file,
)

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
        pyhmmer_results_dict: dict = None,
        custom_pyhmmer_results_dict: dict = None,
        merged_df=None,
        vfdb_results=None,
        card_results=None,
        length_df=None,
        gff_df=None,
        locus_df=None,
        prot_seq_df=None,
        tmrna_flag: bool = False,
        trna_empty: bool = False,
        crispr_count: int = 0,
        coding_table: int = 11,
        mmseqs_flag: bool = True,
        hmm_flag: bool = True,
        custom_hmm_flag: bool = False,
        phanotate_version: str = None,
        pyrodigal_version: str = None,
        pyrodigal_gv_version: str = None,
        pyrodigal_rv_version: str = None,
        trna_version: str = None,
        aragorn_version: str = None,
        minced_version: str = None,
        skip_extra_annotations: bool = False,
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
        merged_df: pl.DataFrame, required
            merged dataframe output
        vfdb_results: pl.DataFrame,
            vfdb dataframe output
        card_results: pl.DataFrame,
            CARD dataframe output
        length_df: pl.DataFrame,
            dataframe with lengths for each input contig
        gff_df: pl.DataFrame,
            dataframe with the output for the gff
        locus_df: pl.DataFrame,
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
        phanotate_version: str or None
            phanotate_version from check_dependencies(); must be set before
            calling process_results() / create_gff()
        pyrodigal_version: str or None
            pyrodigal_version from check_dependencies()
        pyrodigal_gv_version: str or None
            pyrodigal_gv_version from check_dependencies()
        pyrodigal_rv_version: str or None
            pyrodigal_rv_version from check_dependencies()
        trna_version: str or None
            trnascan_version from check_dependencies()
        aragorn_version: str or None
            aragorn_version from check_dependencies()
        minced_version: str or None
            minced_version from check_dependencies()
        prot_seq_df: pl.DataFrame,
            dataframe with protein sequence information for each gene
        skip_extra_annotations: bool
            boolean whether extra annotations are skipped
        reverse_mmseqs2: bool
            use mmseqs2 dbs as target not query
        """
        self.out_dir = out_dir
        self.db_dir = db_dir
        self.prefix = prefix
        self.gene_predictor = gene_predictor
        self.input_fasta = input_fasta
        self.meta_mode = meta_mode
        self.locustag = locustag
        self.pyhmmer_results_dict = pyhmmer_results_dict if pyhmmer_results_dict is not None else {}
        self.custom_pyhmmer_results_dict = custom_pyhmmer_results_dict if custom_pyhmmer_results_dict is not None else {}
        self.merged_df = merged_df if merged_df is not None else pl.DataFrame()
        self.vfdb_results = vfdb_results if vfdb_results is not None else pl.DataFrame()
        self.card_results = card_results if card_results is not None else pl.DataFrame()
        self.length_df = length_df if length_df is not None else pl.DataFrame()
        self.gff_df = gff_df if gff_df is not None else pl.DataFrame()
        self.locus_df = locus_df if locus_df is not None else pl.DataFrame()
        self.prot_seq_df = prot_seq_df if prot_seq_df is not None else pl.DataFrame()
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
        self.pyrodigal_rv_version = pyrodigal_rv_version
        self.trna_version = trna_version
        self.aragorn_version = aragorn_version
        self.minced_version = minced_version
        self.skip_extra_annotations = skip_extra_annotations
        self.reverse_mmseqs2 = reverse_mmseqs2
        # Lazy cache of input FASTA records keyed by contig id.  Built on first
        # access via _get_input_records(); avoids repeatedly re-parsing the
        # FASTA from disk in the downstream writers.
        self._input_records_cache: dict | None = None

    def _get_input_records(self) -> dict:
        """Return ``{contig_id: SeqRecord}`` for the input FASTA, cached.

        Multiple output writers (``create_gff``, ``create_gff_singles``,
        ``split_fasta_singles``) all need to walk the input FASTA.  Parsing
        once and caching avoids redundant disk I/O — and turns
        ``create_gff_singles``' previously O(n × m) re-parse-per-contig into
        an O(n + m) dict lookup.
        """
        if self._input_records_cache is None:
            with open(self.input_fasta) as fa_handle:
                self._input_records_cache = {
                    rec.id: rec for rec in SeqIO.parse(fa_handle, "fasta")
                }
        return self._input_records_cache

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

        if self.gene_predictor in ["prodigal", "prodigal-gv", "pyrodigal-rv"]:
            col_list = ["start", "stop", "strand", "contig", "score", "partial", "gene"]
        else:
            col_list = ["start", "stop", "strand", "contig", "score", "gene"]

        # Read only columns present in file
        cds_df = pl.read_csv(
            cds_file,
            separator="\t",
            skip_rows=1,
            has_header=False,
            infer_schema=False,
        )
        # Assign column names based on how many columns were actually read
        actual_cols = col_list[: cds_df.width]
        cds_df.columns = actual_cols
        cds_df = cds_df.with_columns(pl.col("contig").cast(pl.Utf8))

        ###########################################
        # add the sequence to the df for the genbank conversion later on
        fasta_input_aas_tmp = os.path.join(
            self.out_dir, f"{self.gene_predictor}_aas_tmp.fasta"
        )
        prot_dict = SeqIO.to_dict(SeqIO.parse(fasta_input_aas_tmp, "fasta"))

        # make a copy of cds_df
        self.prot_seq_df = cds_df.clone()

        # to match the output for gff — split gene on first space, keep only the first part
        self.prot_seq_df = self.prot_seq_df.with_columns(
            pl.col("gene").str.splitn(" ", 2).struct.field("field_0").alias("gene")
        )

        # get sequences for each gene in df
        genes_list = self.prot_seq_df["gene"].to_list()
        sequences = []
        for gene in genes_list:
            if gene in prot_dict:
                seq = prot_dict[gene].seq
                sequences.append(str(seq) if seq is not None else "")
            else:
                sequences.append("MA")
        self.prot_seq_df = self.prot_seq_df.with_columns(
            pl.Series("sequence", sequences)
        )

        ##########################################
        # create the tophits_df and write it to file
        if self.mmseqs_flag is True:
            tophits_df = create_mmseqs_tophits(self.out_dir, self.reverse_mmseqs2)

        else:
            # create tophits df
            data = {
                "mmseqs_phrog": ["No_MMseqs"] * cds_df.height,
                "gene": cds_df["gene"].to_list(),
                "mmseqs_alnScore": ["No_MMseqs"] * cds_df.height,
                "mmseqs_seqIdentity": ["No_MMseqs"] * cds_df.height,
                "mmseqs_eVal": ["No_MMseqs"] * cds_df.height,
            }
            tophits_df = pl.DataFrame(data)

        # convert the gene to string for the merge
        cds_df = cds_df.with_columns(pl.col("gene").cast(pl.Utf8))
        tophits_df = tophits_df.with_columns(pl.col("gene").cast(pl.Utf8))
        cds_df = cds_df.filter(pl.col("start").is_not_null())
        cds_df = cds_df.drop_nulls()

        # merge top hits into the cds df
        merged_df = cds_df.join(tophits_df, on="gene", how="left")

        # Adds pyhmmer results if true
        if self.hmm_flag is True:
            merged_df = process_pyhmmer_results(merged_df, self.pyhmmer_results_dict)
        else:
            merged_df = merged_df.with_columns([
                pl.lit("No_PHROGs_HMM").alias("pyhmmer_phrog"),
                pl.lit("No_PHROGs_HMM").alias("pyhmmer_bitscore"),
                pl.lit("No_PHROGs_HMM").alias("pyhmmer_evalue"),
            ])

        # add custom db results
        if self.custom_hmm_flag is True:
            merged_df = process_custom_pyhmmer_results(
                merged_df, self.custom_pyhmmer_results_dict
            )
        else:
            merged_df = merged_df.with_columns([
                pl.lit("No_custom_HMM").alias("custom_hmm_id"),
                pl.lit("No_custom_HMM").alias("custom_hmm_bitscore"),
                pl.lit("No_custom_HMM").alias("custom_hmm_evalue"),
            ])

        # strip off phrog_ for both
        merged_df = merged_df.with_columns([
            pl.col("mmseqs_phrog").str.replace("phrog_", ""),
            pl.col("pyhmmer_phrog").str.replace("phrog_", ""),
        ])

        ############
        # Build the overall "phrog" column.  In mmseqs mode: take mmseqs_phrog,
        # falling back to pyhmmer_phrog for CDS with no mmseqs hit.  Otherwise
        # just copy pyhmmer_phrog.
        if self.mmseqs_flag is True:
            merged_df = merged_df.with_columns(
                pl.when(
                    pl.col("mmseqs_phrog").is_null() & (pl.col("pyhmmer_phrog") != "No_PHROGs_HMM")
                ).then(pl.col("pyhmmer_phrog"))
                .otherwise(pl.col("mmseqs_phrog"))
                .alias("phrog")
            )
        else:  # only need to worry about pyhmmer
            merged_df = merged_df.with_columns(pl.col("pyhmmer_phrog").alias("phrog"))

        # read in phrog annotation file
        phrog_annot_df = pl.read_csv(
            os.path.join(self.db_dir, "phrog_annot_v4.tsv"), separator="\t", infer_schema=False
        )
        phrog_annot_df = phrog_annot_df.with_columns(pl.col("phrog").cast(pl.Utf8))

        # add columns
        if self.gene_predictor == "phanotate":
            method_val = f"ab initio prediction:PHANOTATE:{self.phanotate_version}"
        elif self.gene_predictor == "prodigal":
            method_val = f"ab initio prediction:Pyrodigal:{self.pyrodigal_version}"
        elif self.gene_predictor == "genbank":
            method_val = "CUSTOM"
        elif self.gene_predictor == "prodigal-gv":
            method_val = f"ab initio prediction:Pyrodigal-gv:{self.pyrodigal_gv_version}"
        elif self.gene_predictor == "pyrodigal-rv":
            method_val = f"ab initio prediction:Pyrodigal-rv:{self.pyrodigal_rv_version}"
        else:
            method_val = "UNKNOWN"
        merged_df = merged_df.with_columns([
            pl.lit(method_val).alias("Method"),
            pl.lit("CDS").alias("Region"),
        ])

        # cast mmseqs columns to string
        cols_to_force_string = [
            "mmseqs_phrog",
            "mmseqs_alnScore",
            "mmseqs_seqIdentity",
            "mmseqs_eVal",
        ]

        # Cast mmseqs cols to string and fill nulls with "No_PHROG" in one pass.
        # mmseqs_alnScore is the raw integer string from the mmseqs output (e.g. "82").
        # Differs from pharokka v1.9.1 (pandas): when null rows exist (no-hit CDS),
        # pandas promoted the int64 column to float64 due to NaN, writing "82.0"
        # instead of "82". New polars behaviour keeps the raw integer string.
        # Example: "82.0" (old, when ≥1 CDS had no hit) → "82" (new).
        #
        # The "phrog" column also gets cast + filled here.  "No_PHROG" exists in
        # phrog_annot_v4.tsv with annot="NA"; the explicit "NA"→"hypothetical
        # protein" replacement below handles that case.
        merged_df = merged_df.with_columns([
            pl.col(c).cast(pl.Utf8).fill_null("No_PHROG")
            for c in cols_to_force_string + ["phrog"]
        ])

        # merge phrog annotation (single join) and normalise null/"NA" values
        # in one pass.  phrog_annot_v4.tsv stores "NA" as a literal string for
        # unannotated phrogs — pandas treated "NA" as NaN so fillna() worked,
        # but polars with infer_schema=False keeps it as a string, so we
        # explicitly map both null and "NA" to the user-facing default value.
        merged_df = merged_df.join(phrog_annot_df, on="phrog", how="left").with_columns([
            pl.when(pl.col("annot").is_null() | (pl.col("annot") == "NA"))
              .then(pl.lit("hypothetical protein"))
              .otherwise(pl.col("annot"))
              .alias("annot"),
            pl.when(pl.col("category").is_null() | (pl.col("category") == "NA"))
              .then(pl.lit("unknown function"))
              .otherwise(pl.col("category"))
              .alias("category"),
            pl.col("color").fill_null("None"),
        ])

        # process vfdb results
        (merged_df, vfdb_results) = process_vfdb_results(
            self.out_dir, merged_df, proteins_flag=False, reverse_mmseqs2=self.reverse_mmseqs2
        )
        # process CARD results
        (merged_df, card_results) = process_card_results(
            self.out_dir, merged_df, self.db_dir, proteins_flag=False, reverse_mmseqs2=self.reverse_mmseqs2
        )

        self.merged_df = merged_df
        self.vfdb_results = vfdb_results
        self.card_results = card_results

    def get_contig_name_lengths(self):
        """
        Gets contig name and length in the input fasta file and calculates gc.
        Also adds translation table
        :param fasta_input: input fasta file
        :return: length_df a polars dataframe (to class)
        """

        if self.gene_predictor == "prodigal-gv" or self.gene_predictor == "pyrodigal-rv":
            # define col list
            col_list = [
                "contig",
                "Method",
                "Region",
                "start",
                "stop",
                "score",
                "strand",
                "frame",
                "attributes",
            ]
            # read gff (no fasta output)
            pyrodigal_gv_rv_gff = pl.read_csv(
                os.path.join(self.out_dir, f"{self.gene_predictor}_out.gff"),
                separator="\t",
                has_header=False,
                new_columns=col_list,
                comment_prefix="#",
                truncate_ragged_lines=True,
                infer_schema=False,
            )
            # Split attributes to extract transl_table
            pyrodigal_gv_rv_gff = pyrodigal_gv_rv_gff.with_columns(
                pl.col("attributes").str.splitn("transl_table=", 2).alias("_split1")
            ).with_columns([
                pl.col("_split1").struct.field("field_0").alias("attributes"),
                pl.col("_split1").struct.field("field_1").alias("transl_table"),
            ]).drop("_split1")

            pyrodigal_gv_rv_gff = pyrodigal_gv_rv_gff.with_columns(
                pl.col("transl_table").str.splitn(";conf", 2).alias("_split2")
            ).with_columns([
                pl.col("_split2").struct.field("field_0").alias("transl_table"),
                pl.col("_split2").struct.field("field_1").alias("rest"),
            ]).drop("_split2")

            # drop and then remove duplicates
            pyrodigal_gv_rv_gff = pyrodigal_gv_rv_gff.drop([
                "rest", "Method", "Region", "start", "stop", "score", "strand", "frame", "attributes",
            ])
            pyrodigal_gv_rv_gff = pyrodigal_gv_rv_gff.unique()
            # Convert to a dictionary
            transl_table_dict = dict(
                zip(pyrodigal_gv_rv_gff["contig"].to_list(), pyrodigal_gv_rv_gff["transl_table"].to_list())
            )

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

        with open(self.input_fasta) as fa_handle:
            for record in SeqIO.parse(fa_handle, "fasta"):
                contig_names.append(record.id)
                lengths.append(len(record.seq))
                gc.append(round(gc_fraction(record.seq), 2))
                # pyrodigal-gv lookup from the dict
                if self.gene_predictor == "prodigal-gv" or self.gene_predictor == "pyrodigal-rv":
                    try:
                        transl_table = transl_table_dict[record.id]
                    except Exception:
                        transl_table = "No_CDS_called"

                transl_tables.append(transl_table)

        length_df = pl.DataFrame(
            {
                "contig": contig_names,
                "length": lengths,
                "gc_perc": gc,
                "transl_table": [str(t) for t in transl_tables],
            }
        )
        self.length_df = length_df

    def parse_aragorn(self):
        """
        Parses aragorn output file
        :param out_dir: output directory path
        :param prefix: prefix
        :param length_df: a polars df as output from get_contig_name_lengths()
        :return:- to class - tmrna_flag Boolean whether there is tmrna or not
        """
        aragorn_file = os.path.join(self.out_dir, self.prefix + "_aragorn.txt")
        with open(aragorn_file) as f:
            lines = f.readlines()
        contig_count = self.length_df.height
        tmrna_flag = False
        contig_names = []
        methods = []
        regions = []
        starts = []
        stops = []
        scores = []
        strands = []
        frames = []
        attributes = []
        contig_list = self.length_df["contig"].to_list()
        # if there is only one contig
        if contig_count == 1:
            # if no trnas
            if int(lines[1][0]) == 0:
                tmrna_df = pl.DataFrame(
                    {
                        "contig": [""],
                        "Method": [""],
                        "Region": [""],
                        "start": [""],
                        "stop": [""],
                        "score": [""],
                        "strand": [""],
                        "frame": [""],
                        "attributes": [""],
                    }
                )
            else:
                tmrna_flag = True
                # get all lines with tmrnas
                tmrna_lines = lines[2:]
                for line in tmrna_lines:
                    split = line.split()
                    start_stops = split[2].replace("[", "").replace("]", "").split(",")
                    contig = contig_list[0]
                    method = f"profile:Aragorn:{self.aragorn_version}"
                    region = "tmRNA"
                    start = start_stops[0].replace("c", "")
                    stop = start_stops[1]
                    score = "."
                    strand = "."
                    frame = "."
                    tag_peptide = split[3].replace(",", "..")
                    tag_peptide_seq = split[4]
                    attribute = (
                        "product=transfer-messenger RNA SsrA;tag_peptide="
                        + tag_peptide
                        + ";note=tag peptide sequence: "
                        + tag_peptide_seq
                    )
                    contig_names.append(contig)
                    methods.append(method)
                    regions.append(region)
                    starts.append(start)
                    stops.append(stop)
                    scores.append(score)
                    strands.append(strand)
                    frames.append(frame)
                    attributes.append(attribute)
                tmrna_df = pl.DataFrame(
                    {
                        "contig": contig_names,
                        "Method": methods,
                        "Region": regions,
                        "start": starts,
                        "stop": stops,
                        "score": scores,
                        "strand": strands,
                        "frame": frames,
                        "attributes": attributes,
                    }
                )
        # two or more contigs
        else:
            i = 0  # line counter
            j = 0  # contig counter
            for line in lines:
                if (
                    line[0] == ">" and "sensitivity" not in line
                ):
                    if "0 genes found" not in lines[i + 1]:
                        tmrna_flag = True
                        tmrna_count = int(lines[i + 1][0])
                        for k in range(tmrna_count):
                            tmrna_line = lines[i + 2 + k]
                            split = tmrna_line.split()
                            start_stops = (
                                split[2].replace("[", "").replace("]", "").split(",")
                            )
                            contig = contig_list[j]
                            method = f"profile:Aragorn:{self.aragorn_version}"
                            region = "tmRNA"
                            start = start_stops[0].replace("c", "")
                            stop = start_stops[1]
                            score = "."
                            strand = "."
                            frame = "."
                            tag_peptide = split[3].replace(",", "..")
                            tag_peptide_seq = split[4]
                            attribute = (
                                "product=transfer-messenger RNA SsrA;tag_peptide="
                                + tag_peptide
                                + ";note=tag peptide sequence: "
                                + tag_peptide_seq
                            )
                            contig_names.append(contig)
                            methods.append(method)
                            regions.append(region)
                            starts.append(start)
                            stops.append(stop)
                            scores.append(score)
                            strands.append(strand)
                            frames.append(frame)
                            attributes.append(attribute)
                    j += 1
                i += 1
            tmrna_df = pl.DataFrame(
                {
                    "contig": contig_names,
                    "Method": methods,
                    "Region": regions,
                    "start": starts,
                    "stop": stops,
                    "score": scores,
                    "strand": strands,
                    "frame": frames,
                    "attributes": attributes,
                }
            )
        tmrna_df.write_csv(
            os.path.join(self.out_dir, self.prefix + "_aragorn.gff"),
            separator="\t",
            include_header=False,
            quote_style="never",
        )
        self.tmrna_flag = tmrna_flag

    def create_gff(self):
        """
        Creates the pharokka.gff file
        """

        # create locus tag
        if self.locustag == "Random":
            self.locustag = "".join(
                random.choice(string.ascii_uppercase) for _ in range(8)
            )

        # get all contigs
        contigs = self.length_df["contig"].cast(pl.Utf8)
        contigs_list = contigs.to_list()
        self.merged_df = self.merged_df.with_columns(pl.col("contig").cast(pl.Utf8))

        # add the translation table
        transl_table_df = self.length_df.drop(["length", "gc_perc"])
        transl_table_df = transl_table_df.with_columns(pl.col("contig").cast(pl.Utf8))
        self.merged_df = self.merged_df.join(transl_table_df, how="left", on="contig")

        ############ locus tag #########
        if self.meta_mode:
            subset_dfs = []
            for contig in contigs_list:
                subset_df = self.merged_df.filter(pl.col("contig") == contig)
                subset_df = subset_df.with_row_index("_idx").with_columns(
                    (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("count")
                ).drop("_idx")
                subset_dfs.append(subset_df)
            self.merged_df = pl.concat(subset_dfs)
        else:
            self.merged_df = self.merged_df.with_row_index("_idx").with_columns(
                (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("count")
            ).drop("_idx")

        # get the locus tag
        if not self.meta_mode:
            self.merged_df = self.merged_df.with_columns(
                (pl.lit(self.locustag + "_CDS_") + pl.col("count")).alias("locus_tag")
            )
        else:
            self.merged_df = self.merged_df.with_columns(
                (pl.col("contig") + pl.lit("_CDS_") + pl.col("count")).alias("locus_tag")
            )

        # Propagate the locus_tag to prot_seq_df so that convert_gff_to_gbk
        # can key protein translations by locus_tag rather than by row
        # position (which silently produced wrong /translation values if the
        # prot_seq_df rows ever diverged from GFF feature order).
        #
        # prot_seq_df.gene was split on first space upstream (in process_results)
        # whereas merged_df.gene retains the full pre-split string, so we
        # normalise merged_df's gene with the same splitn before joining.
        self.prot_seq_df = self.prot_seq_df.with_columns(pl.col("gene").cast(pl.Utf8))
        _locus_lookup = self.merged_df.select([
            pl.col("gene").cast(pl.Utf8).str.splitn(" ", 2).struct.field("field_0").alias("gene"),
            "locus_tag",
        ])
        self.prot_seq_df = self.prot_seq_df.join(_locus_lookup, on="gene", how="left")

        #########
        # Rearrange start and stop so that start is always less than stop for GFF.
        # GFF3 requires start ≤ end (purely positional); strand carries direction.
        #########
        self.merged_df = self.merged_df.with_columns([
            pl.when(pl.col("strand") == "-").then(pl.col("stop")).otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("strand") == "-").then(pl.col("start")).otherwise(pl.col("stop")).alias("stop"),
        ])

        # Capture locus_df AFTER the GFF coordinate swap so that start ≤ stop in
        # both single-phage and meta mode.  locus_df is used by
        # write_tophits_vfdb_card() to supply start/stop/strand for
        # top_hits_vfdb.tsv and top_hits_card.tsv.
        #
        # v1.9.1 bug (now fixed): in meta mode, `locus_df = pd.concat(subset_dfs)`
        # created a brand-new pandas DataFrame, breaking the reference to
        # self.merged_df before the GFF swap ran.  Consequently locus_df kept
        # biological-order coordinates (start > stop for negative-strand CDS), so
        # top_hits_vfdb/card in meta mode incorrectly reported e.g.
        # start=40487, stop=40038 for a neg-strand hit instead of the correct
        # start=40038, stop=40487.  Single-phage mode was unaffected because
        # `locus_df = self.merged_df` was a reference alias that shared the
        # in-place mutation.  Both modes now capture after the swap.
        self.locus_df = self.merged_df

        # set frame to be 0
        self.merged_df = self.merged_df.with_columns(pl.lit("0").alias("frame"))

        # create attributes
        self.merged_df = self.merged_df.with_columns(
            (
                pl.lit("ID=") + pl.col("locus_tag").cast(pl.Utf8)
                + pl.lit(";transl_table=") + pl.col("transl_table").cast(pl.Utf8)
                + pl.lit(";phrog=") + pl.col("phrog").fill_null("No_PHROG").cast(pl.Utf8)
                + pl.lit(";locus_tag=") + pl.col("locus_tag").cast(pl.Utf8)
                + pl.lit(";function=") + pl.col("category").cast(pl.Utf8)
                + pl.lit(";product=") + pl.col("annot").cast(pl.Utf8)
            ).alias("attributes")
        )

        if "partial" in self.merged_df.columns:
            self.merged_df = self.merged_df.with_columns(
                (pl.col("attributes") + pl.lit(";partial=") + pl.col("partial").cast(pl.Utf8)).alias("attributes")
            )

        # adds custom hmm database annotations
        self.merged_df = self.merged_df.with_columns(
            pl.when(pl.col("custom_hmm_id") != "No_custom_HMM")
            .then(
                pl.col("attributes").cast(pl.Utf8)
                + pl.lit(";custom_annotation=") + pl.col("custom_hmm_id").cast(pl.Utf8)
            )
            .otherwise(pl.col("attributes"))
            .alias("attributes")
        )
        # adds VFDB
        self.merged_df = self.merged_df.with_columns(
            pl.when(pl.col("vfdb_short_name") != "None")
            .then(
                pl.col("attributes").cast(pl.Utf8)
                + pl.lit(";vfdb_short_name=") + pl.col("vfdb_short_name").cast(pl.Utf8)
                + pl.lit(";vfdb_description=") + pl.col("vfdb_description").cast(pl.Utf8)
                + pl.lit(";vfdb_species=") + pl.col("vfdb_species").cast(pl.Utf8)
            )
            .otherwise(pl.col("attributes"))
            .alias("attributes")
        )
        # adds CARD
        self.merged_df = self.merged_df.with_columns(
            pl.when(pl.col("CARD_short_name") != "None")
            .then(
                pl.col("attributes").cast(pl.Utf8)
                + pl.lit(";CARD_short_name=") + pl.col("CARD_short_name").cast(pl.Utf8)
                + pl.lit(";AMR_Gene_Family=") + pl.col("AMR_Gene_Family").cast(pl.Utf8)
                + pl.lit(";CARD_species=") + pl.col("CARD_species").cast(pl.Utf8)
            )
            .otherwise(pl.col("attributes"))
            .alias("attributes")
        )

        # get gff dataframe in correct order
        gff_df = self.merged_df.select([
            "contig", "Method", "Region", "start", "stop", "score", "strand", "frame", "attributes",
        ])

        # change start and stop to int
        gff_df = gff_df.with_columns([
            pl.col("start").cast(pl.Int64),
            pl.col("stop").cast(pl.Int64),
        ])

        ### trnas
        if self.skip_extra_annotations is False:
            col_list = [
                "contig", "Method", "Region", "start", "stop", "score", "strand", "frame", "attributes",
            ]
            trna_empty = is_trna_empty(self.out_dir)
            if not trna_empty:
                trna_df = pl.read_csv(
                    os.path.join(self.out_dir, "trnascan_out.gff"),
                    separator="\t",
                    has_header=False,
                    new_columns=col_list,
                    comment_prefix="#",
                    truncate_ragged_lines=True,
                    infer_schema=False,
                )
                trna_df = trna_df.with_columns(pl.col("contig").cast(pl.Utf8))
                trna_df = trna_df.with_columns(
                    pl.lit(f"profile:tRNAscan-SE:{self.trna_version}").alias("Method")
                )

                # Extract the tRNAscan ID (e.g. "MW460250_1.trna3") from
                # attributes — used to key the anticodon-positions lookup
                # below.  Raw attribute string starts with "ID=<id>;..."
                trna_df = trna_df.with_columns(
                    pl.col("attributes")
                      .str.extract(r"ID=([^;]+)", 1)
                      .alias("_trnascan_id")
                )

                # index hack if meta mode
                if self.meta_mode:
                    subset_dfs = []
                    for contig in contigs_list:
                        subset_df = trna_df.filter(pl.col("contig") == contig)
                        subset_df = subset_df.filter(
                            (pl.col("Region") == "tRNA") | (pl.col("Region") == "pseudogene")
                        )
                        subset_df = subset_df.with_row_index("_idx").with_columns(
                            (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("_count")
                        ).drop("_idx")
                        subset_df = subset_df.with_columns(
                            (pl.lit(contig + "_tRNA_") + pl.col("_count")).alias("locus_tag")
                        ).drop("_count")
                        subset_dfs.append(subset_df)
                    trna_df = pl.concat(subset_dfs)
                else:
                    trna_df = trna_df.filter(
                        (pl.col("Region") == "tRNA") | (pl.col("Region") == "pseudogene")
                    )
                    trna_df = trna_df.with_row_index("_idx").with_columns(
                        (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("_count")
                    ).drop("_idx")
                    trna_df = trna_df.with_columns(
                        (pl.lit(self.locustag + "_tRNA_") + pl.col("_count")).alias("locus_tag")
                    ).drop("_count")

                trna_df = trna_df.with_columns([
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                ])

                # Split attributes to extract isotype and anticodon
                trna_df = trna_df.with_columns(
                    pl.col("attributes").str.splitn(";isotype=", 2).alias("_s1")
                ).with_columns([
                    pl.col("_s1").struct.field("field_0").alias("attributes"),
                    pl.col("_s1").struct.field("field_1").alias("isotypes"),
                ]).drop("_s1")

                trna_df = trna_df.with_columns(
                    pl.col("isotypes").str.splitn(";anticodon=", 2).alias("_s2")
                ).with_columns([
                    pl.col("_s2").struct.field("field_0").alias("isotypes"),
                    pl.col("_s2").struct.field("field_1").alias("anticodon"),
                ]).drop("_s2")

                trna_df = trna_df.with_columns(
                    pl.col("anticodon").str.splitn(";gene_biotype", 2).alias("_s3")
                ).with_columns([
                    pl.col("_s3").struct.field("field_0").alias("anticodon"),
                    pl.col("_s3").struct.field("field_1").alias("rest"),
                ]).drop("_s3")

                # Extract anticodon positions once, up-front.
                anticodon_positions = extract_anticodon_positions(self.out_dir)

                if "note" not in trna_df.columns:
                    trna_df = trna_df.with_columns(pl.lit("").alias("note"))

                # ─── Vectorised tRNA isotype/codon/note/anticodon_gb logic ───
                # Old code: per-row Python loop building isotypes/note/codon
                # /anticodon_gb lists, then re-attaching them as polars Series.
                # New code: native polars expressions in a single with_columns.
                #
                # Per-row semantics (per the original loop):
                #   isotype == "Undet"  → isotypes := "OTHER", note := "Undetermined tRNA",
                #                         codon := "", anticodon_gb := None
                #   isotype == "Sup"    → isotypes := "TERM",  note := "Suppressor tRNA",
                #                         codon := get_codon_from_anticodon(anticodon),
                #                         anticodon_gb := "(pos:S..E,aa:TERM,seq:ANT)"
                #   else                → isotypes/note unchanged,
                #                         codon := get_codon_from_anticodon(anticodon),
                #                         anticodon_gb := "(pos:S..E,aa:ISO,seq:ANT)"
                # where (S, E) comes from extract_anticodon_positions() if the
                # row's tRNA ID is present, else the row's own (start, stop).

                # 1) Build a positions DataFrame keyed by the tRNAscan ID.
                # This replaces the v1.9.1 row-position-based lookup, which
                # silently attached anticodon positions to the wrong tRNA when
                # the GFF and .sec files were sorted differently (common — the
                # GFF is in genomic order; the .sec is in tRNA-ID order).
                if anticodon_positions:
                    positions_df = pl.DataFrame({
                        "_trnascan_id": list(anticodon_positions.keys()),
                        "_pos_start":   [s for s, _ in anticodon_positions.values()],
                        "_pos_end":     [e for _, e in anticodon_positions.values()],
                    }).with_columns([
                        pl.col("_pos_start").cast(pl.Int64),
                        pl.col("_pos_end").cast(pl.Int64),
                    ])
                    trna_df = trna_df.join(positions_df, on="_trnascan_id", how="left")
                else:
                    # No .sec file / no positions parsed — fall back to start/stop.
                    trna_df = trna_df.with_columns([
                        pl.lit(None, dtype=pl.Int64).alias("_pos_start"),
                        pl.lit(None, dtype=pl.Int64).alias("_pos_end"),
                    ])

                # 2) Compute the codon column using Biopython's reverse-complement
                # + transcribe.  ``map_elements`` evaluates per-row; the safe
                # wrapper protects against the Undet case where anticodon may be
                # malformed (Biopython would raise) — when/otherwise evaluates
                # both branches in polars so we cannot rely on the Undet mask
                # alone to skip the call.
                def _safe_codon(s):
                    if not isinstance(s, str) or not s.strip():
                        return ""
                    try:
                        return get_codon_from_anticodon(s)
                    except (ValueError, TypeError):
                        return ""

                # 3) Build all derived columns in a single pass.  All when() /
                # otherwise() expressions read the *original* isotypes value
                # since polars evaluates with_columns inputs in parallel against
                # the current frame state.
                _isotypes_remapped = (
                    pl.when(pl.col("isotypes") == "Undet").then(pl.lit("OTHER"))
                      .when(pl.col("isotypes") == "Sup").then(pl.lit("TERM"))
                      .otherwise(pl.col("isotypes"))
                )

                trna_df = trna_df.with_columns([
                    # note: assigned only for Undet / Sup, unchanged otherwise
                    pl.when(pl.col("isotypes") == "Undet").then(pl.lit("Undetermined tRNA"))
                      .when(pl.col("isotypes") == "Sup").then(pl.lit("Suppressor tRNA"))
                      .otherwise(pl.col("note"))
                      .alias("note"),
                    # codon: empty for Undet, computed via Biopython otherwise
                    pl.when(pl.col("isotypes") == "Undet")
                      .then(pl.lit(""))
                      .otherwise(
                          pl.col("anticodon").map_elements(_safe_codon, return_dtype=pl.Utf8)
                      )
                      .alias("codon"),
                    # anticodon_gb: None for Undet; otherwise pos:start..end,aa:ISO,seq:ANT
                    pl.when(pl.col("isotypes") == "Undet")
                      .then(pl.lit(None, dtype=pl.Utf8))
                      .otherwise(
                          pl.format(
                              "(pos:{}..{},aa:{},seq:{})",
                              pl.coalesce([pl.col("_pos_start"), pl.col("start")]),
                              pl.coalesce([pl.col("_pos_end"),   pl.col("stop")]),
                              _isotypes_remapped,
                              pl.col("anticodon"),
                          )
                      )
                      .alias("anticodon_gb"),
                    # isotypes: remapped Undet→OTHER, Sup→TERM, else unchanged
                    _isotypes_remapped.alias("isotypes"),
                ])

                # Drop temp columns from the positions join + ID extraction.
                trna_df = trna_df.drop(["_pos_start", "_pos_end", "_trnascan_id"])

                trna_df = trna_df.with_columns([
                    (
                        pl.lit("tRNA-") + pl.col("isotypes")
                        + pl.when(pl.col("codon") != "")
                        .then(pl.lit("(") + pl.col("codon") + pl.lit(")"))
                        .otherwise(pl.lit(""))
                    ).alias("trna_gene"),
                    (
                        pl.lit("tRNA-") + pl.col("isotypes")
                        + pl.when(pl.col("codon") != "")
                        .then(pl.lit("(") + pl.col("codon") + pl.lit(")"))
                        .otherwise(pl.lit(""))
                    ).alias("trna_product"),
                ])

                if "anticodon_gb" not in trna_df.columns:
                    trna_df = trna_df.with_columns(pl.lit(None).alias("anticodon_gb"))

                trna_df = trna_df.drop(["attributes"])
                trna_df = trna_df.with_columns(
                    (
                        pl.lit("ID=") + pl.col("locus_tag")
                        + pl.lit(";locus_tag=") + pl.col("locus_tag")
                        + pl.lit(";gene=") + pl.col("trna_gene").cast(pl.Utf8)
                        + pl.lit(";product=") + pl.col("trna_product").cast(pl.Utf8)
                        + pl.when(
                            pl.col("anticodon_gb").is_not_null()
                            & (pl.col("anticodon_gb").cast(pl.Utf8) != pl.lit(""))
                        ).then(pl.lit(";anticodon=") + pl.col("anticodon_gb").cast(pl.Utf8))
                        .otherwise(pl.lit(""))
                        + pl.when(
                            pl.col("note").is_not_null()
                            & (pl.col("note").cast(pl.Utf8) != pl.lit(""))
                        ).then(pl.lit(";note=") + pl.col("note").cast(pl.Utf8))
                        .otherwise(pl.lit(""))
                    ).alias("attributes")
                )

                # Keep "codon" — old pharokka's pandas to_csv wrote it as a 10th
                # column in the GFF (e.g. "\tAUG" for tRNA-Met). Non-tRNA rows
                # get null which writes as an empty field (trailing tab).
                # Explicitly select in the correct column order so "codon" comes
                # after "attributes" when written to GFF.
                trna_df = trna_df.select([
                    "contig", "Method", "Region", "start", "stop",
                    "score", "strand", "frame", "attributes", "codon",
                ])

            ### crisprs
            crispr_count = get_crispr_count(self.out_dir, self.prefix)
            # add to gff if > 0
            if crispr_count > 0:
                minced_df = pl.read_csv(
                    os.path.join(self.out_dir, self.prefix + "_minced.gff"),
                    separator="\t",
                    has_header=False,
                    new_columns=col_list,
                    comment_prefix="#",
                    infer_schema=False,
                )
                minced_df = minced_df.with_columns(pl.col("contig").cast(pl.Utf8))
                minced_df = minced_df.with_columns(
                    pl.lit(f"nucleotide motif:MinCED:{self.minced_version}").alias("Method")
                )
                minced_df = minced_df.with_columns([
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                ])

                # split attributes
                minced_df = minced_df.with_columns(
                    pl.col("attributes").str.splitn(";rpt_unit_seq=", 2).alias("_s1")
                ).with_columns([
                    pl.col("_s1").struct.field("field_0").alias("attributes"),
                    pl.col("_s1").struct.field("field_1").alias("rpt_unit_seq"),
                ]).drop("_s1")

                minced_df = minced_df.with_columns(
                    pl.col("attributes").str.splitn(";rpt_family=", 2).alias("_s2")
                ).with_columns([
                    pl.col("_s2").struct.field("field_0").alias("attributes"),
                    pl.col("_s2").struct.field("field_1").alias("rpt_family"),
                ]).drop("_s2")

                minced_df = minced_df.with_columns(
                    pl.col("attributes").str.splitn(";rpt_type=", 2).alias("_s3")
                ).with_columns([
                    pl.col("_s3").struct.field("field_0").alias("attributes"),
                    pl.col("_s3").struct.field("field_1").alias("rpt_type"),
                ]).drop("_s3")

                minced_df = minced_df.drop(["attributes"])

                # index hack if meta mode
                if self.meta_mode:
                    subset_dfs = []
                    for contig in contigs_list:
                        subset_df = minced_df.filter(pl.col("contig") == contig)
                        subset_df = subset_df.with_row_index("_idx").with_columns(
                            (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("_count")
                        ).drop("_idx")
                        subset_df = subset_df.with_columns(
                            (pl.lit(contig + "_CRISPR_") + pl.col("_count").str.zfill(4)).alias("locus_tag")
                        ).drop("_count")
                        subset_dfs.append(subset_df)
                    minced_df = pl.concat(subset_dfs)
                else:
                    minced_df = minced_df.with_row_index("_idx").with_columns(
                        (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("_count")
                    ).drop("_idx")
                    minced_df = minced_df.with_columns(
                        (pl.lit(self.locustag + "_CRISPR_") + pl.col("_count")).alias("locus_tag")
                    ).drop("_count")

                minced_df = minced_df.with_columns(
                    (
                        pl.lit("ID=") + pl.col("locus_tag")
                        + pl.lit(";rpt_type=") + pl.col("rpt_type").cast(pl.Utf8)
                        + pl.lit(";rpt_family=") + pl.col("rpt_family").cast(pl.Utf8)
                        + pl.lit(";rpt_unit_seq=") + pl.col("rpt_unit_seq").cast(pl.Utf8)
                        + pl.lit(";rpt_unit_range=") + pl.col("start").cast(pl.Utf8)
                        + pl.lit("..") + pl.col("stop").cast(pl.Utf8)
                        + pl.lit(";locus_tag=") + pl.col("locus_tag")
                    ).alias("attributes")
                )
                minced_df = minced_df.drop(["rpt_unit_seq", "rpt_family", "rpt_type", "locus_tag"])

            ### tmrna
            if self.tmrna_flag:
                tmrna_df = pl.read_csv(
                    os.path.join(self.out_dir, self.prefix + "_aragorn.gff"),
                    separator="\t",
                    has_header=False,
                    new_columns=col_list,
                    comment_prefix="#",
                    truncate_ragged_lines=True,
                    infer_schema=False,
                )
                tmrna_df = tmrna_df.with_columns(pl.col("contig").cast(pl.Utf8))
                tmrna_df = tmrna_df.with_columns([
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                ])

                # index hack if meta mode
                if self.meta_mode:
                    subset_dfs = []
                    for contig in contigs_list:
                        subset_df = tmrna_df.filter(pl.col("contig") == contig)
                        subset_df = subset_df.with_row_index("_idx").with_columns(
                            (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("_count")
                        ).drop("_idx")
                        subset_df = subset_df.with_columns(
                            (pl.lit(contig + "_tmRNA_") + pl.col("_count")).alias("locus_tag")
                        ).drop("_count")
                        subset_dfs.append(subset_df)
                    tmrna_df = pl.concat(subset_dfs)
                else:
                    tmrna_df = tmrna_df.with_row_index("_idx").with_columns(
                        (pl.col("_idx") + 1).cast(pl.Utf8).str.zfill(4).alias("_count")
                    ).drop("_idx")
                    tmrna_df = tmrna_df.with_columns(
                        (pl.lit(self.locustag + "_tmRNA_") + pl.col("_count")).alias("locus_tag")
                    ).drop("_count")

                tmrna_df = tmrna_df.with_columns(
                    (
                        pl.lit("ID=") + pl.col("locus_tag")
                        + pl.lit(";") + pl.col("attributes").cast(pl.Utf8)
                        + pl.lit(";locus_tag=") + pl.col("locus_tag")
                    ).alias("attributes")
                )
                tmrna_df = tmrna_df.drop(["locus_tag"])

        # write header of final gff files
        with open(os.path.join(self.out_dir, self.prefix + ".gff"), "w") as f:
            f.write("##gff-version 3\n")
            for row in self.length_df.iter_rows(named=True):
                f.write(
                    "##sequence-region "
                    + row["contig"]
                    + " 1 "
                    + str(row["length"])
                    + "\n"
                )

        # combine dfs depending on whether the elements were detected
        if self.skip_extra_annotations is True:
            df_list = [gff_df]
        else:
            if trna_empty is True and self.tmrna_flag is False and crispr_count == 0:
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
            else:  # all detected
                df_list = [gff_df, trna_df, tmrna_df, minced_df]

        total_gff = pl.concat(df_list, how="diagonal")

        # ensure that the start and stops are integer
        total_gff = total_gff.with_columns([
            pl.col("start").cast(pl.Int64),
            pl.col("stop").cast(pl.Int64),
        ])

        # sorts all features by start
        total_gff = total_gff.sort(["contig", "start"])

        # write final gff to file
        # Old pharokka used pandas to_csv on a DataFrame that had a "codon"
        # column (kept from trna_df). For tRNA rows it wrote e.g. "\tAUG" as a
        # 10th field; for non-tRNA rows the null codon produced a trailing "\t".
        # We replicate that by keeping "codon" in trna_df and using diagonal
        # concat so non-tRNA rows get null → empty string → trailing tab.
        with open(os.path.join(self.out_dir, self.prefix + ".gff"), "a") as f:
            total_gff.write_csv(f, separator="\t", include_header=False, quote_style="never", null_value="")

        # write fasta on the end
        with open(os.path.join(self.out_dir, self.prefix + ".gff"), "a") as f:
            f.write("##FASTA\n")
            for record in self._get_input_records().values():
                f.write(f">{record.id}\n")
                sequence = record.seq
                chunk_size = 60
                for i in range(0, len(sequence), chunk_size):
                    f.write(str(sequence[i: i + chunk_size]) + "\n")

        self.gff_df = gff_df
        self.total_gff = total_gff
        if self.skip_extra_annotations is False:
            self.trna_empty = trna_empty
            self.crispr_count = crispr_count
        else:
            self.trna_empty = True
            self.crispr_count = 0
            self.tmrna_flag = False

    def create_tbl(self):
        """
        Creates the pharokka.tbl file
        """

        col_list = [
            "contig", "Method", "Region", "start", "stop", "score", "strand", "frame", "attributes",
        ]

        ### trnas
        if self.trna_empty is False:
            trna_df = self.total_gff.filter(
                pl.col("Method") == f"profile:tRNAscan-SE:{self.trna_version}"
            )
            trna_df = parse_attributes_column(trna_df)
            trna_df = trna_df.with_columns([
                pl.col("contig").cast(pl.Utf8),
                pl.col("start").cast(pl.Int64),
                pl.col("stop").cast(pl.Int64),
            ])
            if "anticodon" not in trna_df.columns:
                trna_df = trna_df.with_columns(pl.lit(None).alias("anticodon"))

        #### CRISPRs
        if self.crispr_count > 0:
            crispr_df = self.total_gff.filter(pl.col("Region") == "repeat_region")
            crispr_df = parse_attributes_column(crispr_df)
            crispr_df = crispr_df.with_columns([
                pl.col("contig").cast(pl.Utf8),
                pl.col("start").cast(pl.Int64),
                pl.col("stop").cast(pl.Int64),
            ])

        ### TMRNAs
        if self.tmrna_flag is True:
            tmrna_df = self.total_gff.filter(pl.col("Region") == "tmRNA")
            tmrna_df = parse_attributes_column(tmrna_df)
            tmrna_df = tmrna_df.with_columns([
                pl.col("contig").cast(pl.Utf8),
                pl.col("start").cast(pl.Int64),
                pl.col("stop").cast(pl.Int64),
            ])

        with open(os.path.join(self.out_dir, self.prefix + ".tbl"), "w") as f:
            for row in self.length_df.iter_rows(named=True):
                contig = str(row["contig"])
                contig_length = int(row["length"])
                f.write(">Feature " + contig + "\n")
                subset_df = self.merged_df.filter(pl.col("contig") == contig)
                # drop transl_table and partial column
                if "transl_table" in subset_df.columns:
                    subset_df = subset_df.drop(["transl_table"])
                if "partial" in subset_df.columns:
                    subset_df = subset_df.drop(["partial"])
                subset_df = parse_attributes_column(subset_df)
                subset_df = subset_df.with_columns([
                    pl.col("start").cast(pl.Int64),
                    pl.col("stop").cast(pl.Int64),
                ])
                for srow in subset_df.iter_rows(named=True):
                    start = str(srow["start"])
                    stop = str(srow["stop"])
                    codon_start = None

                    if srow["strand"] == "-":
                        start = str(srow["stop"])
                        stop = str(srow["start"])

                    if "partial" in srow and srow.get("partial") is not None:
                        if srow["strand"] == "+" and srow["partial"] == "10":
                            start = "<" + str(srow["start"])
                            if srow["start"] > 1:
                                codon_start = str(srow["start"])
                                start = "<1"
                        elif srow["strand"] == "+" and srow["partial"] == "01":
                            stop = ">" + str(srow["stop"])
                        elif srow["strand"] == "-" and srow["partial"] == "01":
                            start = "<" + str(srow["stop"])
                            if contig_length - srow["stop"] > 0:
                                codon_start = (
                                    contig_length - srow["stop"] + 1
                                )
                                start = "<" + str(contig_length)
                                if codon_start > 3:
                                    logger.error(
                                        "Error: codon_start can not be greater than 3. Please raise an issue on GitHub with your genome. \n"
                                    )
                        elif srow["strand"] == "-" and srow["partial"] == "10":
                            stop = ">" + str(srow["start"])

                    f.write(f"{start}\t{stop}\t{srow['Region']}\n")
                    f.write(f"\t\t\tproduct\t{srow['annot']}\n")
                    f.write(f"\t\t\tfunction\t{srow['function']}\n")
                    f.write(f"\t\t\tinference\t{srow['Method']}\n")
                    f.write(f"\t\t\ttransl_table\t{srow['transl_table']}\n")
                    if codon_start is not None:
                        f.write(f"\t\t\tcodon_start\t{codon_start}\n")
                if not self.trna_empty:
                    subset_trna_df = trna_df.filter(pl.col("contig") == contig)
                    for trow in subset_trna_df.iter_rows(named=True):
                        start = str(trow["start"])
                        stop = str(trow["stop"])
                        if trow["strand"] == "-":
                            start = str(trow["stop"])
                            stop = str(trow["start"])
                        f.write(f"{start}\t{stop}\t{trow['Region']}\n")
                        f.write(f"\t\t\tgene\t{trow['gene']}\n")
                        f.write(f"\t\t\tproduct\t{trow['product']}\n")
                        f.write(f"\t\t\tinference\t{trow['Method']}\n")
                        if trow.get("anticodon") is not None:
                            f.write(f"\t\t\tanticodon\t{trow['anticodon']}\n")
                        if trow.get("note") is not None:
                            f.write(f"\t\t\tnote\t{trow['note']}\n")
                if self.crispr_count > 0:
                    subset_crispr_df = crispr_df.filter(pl.col("contig") == contig)
                    for crow in subset_crispr_df.iter_rows(named=True):
                        start = str(crow["start"])
                        stop = str(crow["stop"])
                        if crow["strand"] == "-":
                            start = str(crow["stop"])
                            stop = str(crow["start"])
                        f.write(f"{start}\t{stop}\t{crow['Region']}\n")
                        f.write(f"\t\t\tinference\t{crow['Method']}\n")
                        f.write(f"\t\t\trpt_family\t{crow['rpt_family']}\n")
                        f.write(f"\t\t\trpt_type\t{crow['rpt_type']}\n")
                        f.write(f"\t\t\trpt_unit_range\t{crow['rpt_unit_range']}\n")
                        f.write(f"\t\t\trpt_unit_seq\t{crow['rpt_unit_seq']}\n")
                if self.tmrna_flag:
                    subset_tmrna_df = tmrna_df.filter(pl.col("contig") == contig)
                    for tmrow in subset_tmrna_df.iter_rows(named=True):
                        start = str(tmrow["start"])
                        stop = str(tmrow["stop"])
                        if tmrow["strand"] == "-":
                            start = str(tmrow["stop"])
                            stop = str(tmrow["start"])
                        f.write(f"{start}\t{stop}\ttmRNA\n")
                        f.write(f"\t\t\tinference\t{tmrow['Method']}\n")
                        f.write(f"\t\t\tproduct\t{tmrow['product']}\n")
                        f.write(f"\t\t\ttag_peptide\t{tmrow['tag_peptide']}\n")
                        f.write(f"\t\t\tnote\t{tmrow['note']}\n")

    def create_gff_singles(self):
        """
        Creates the single gff3 files for each input contig in meta
        """

        single_gff_dir = os.path.join(self.out_dir, "single_gffs")
        check_and_create_directory(single_gff_dir)

        for row in self.length_df.iter_rows(named=True):
            contig = row["contig"]
            with open(os.path.join(single_gff_dir, contig + ".gff"), "w") as f:
                f.write("##gff-version 3\n")
                f.write(
                    "##sequence-region " + contig + " 1 " + str(row["length"]) + "\n"
                )

            subset_df = self.total_gff.filter(pl.col("contig") == contig)

            with open(os.path.join(single_gff_dir, contig + ".gff"), "a") as f:
                subset_df.write_csv(f, separator="\t", include_header=False, quote_style="never")

            with open(os.path.join(single_gff_dir, contig + ".gff"), "a") as f:
                f.write("##FASTA\n")
                # O(1) dict lookup — previously this re-parsed the entire FASTA
                # once per contig (O(n × m) in total).
                dna_record = self._get_input_records().get(contig)
                if dna_record is not None:
                    SeqIO.write(dna_record, f, "fasta")

    def convert_singles_gff_to_gbk(self):
        """
        Converts all single gffs to gbks
        """

        single_gff_dir = os.path.join(self.out_dir, "single_gffs")
        single_gbk_dir = os.path.join(self.out_dir, "single_gbks")
        check_and_create_directory(single_gbk_dir)

        split_fasta_dir = os.path.join(self.out_dir, "input_split_tmp")

        for i, row in enumerate(self.length_df.iter_rows(named=True)):
            fasta_file = os.path.join(
                split_fasta_dir, "input_subprocess" + str(i + 1) + ".fasta"
            )
            contig = row["contig"]
            convert_gff_to_gbk(
                fasta_file, single_gff_dir, single_gbk_dir, contig, self.prot_seq_df
            )

    def split_fasta_singles(self):
        """Splits the input fasta into separate single fasta files for output based on contig names"""

        single_fastas = os.path.join(self.out_dir, "single_fastas")
        check_and_create_directory(single_fastas)
        for dna_record in self._get_input_records().values():
            contig = dna_record.id
            with open(os.path.join(single_fastas, contig + ".fasta"), "w") as f:
                SeqIO.write(dna_record, f, "fasta")

    def split_faas_singles(self):
        """Splits the .faa fasta into separate single fasta files for output based on contig names"""

        single_faas = os.path.join(self.out_dir, "single_faas")
        check_and_create_directory(single_faas)
        faa_file = os.path.join(self.out_dir, f"{self.gene_predictor}.faa")
        with open(faa_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                protein_id = record.id[:-9]
                with open(os.path.join(single_faas, f"{protein_id}.faa"), "a") as f:
                    SeqIO.write(record, f, "fasta")

    def write_tophits_vfdb_card(self):
        """
        Outputs top_hits_vfdb.tsv and top_hits_card.tsv
        """
        ######################################
        ##### update vfdb with locus tag #####
        ######################################
        # locus_df is self.merged_df (set in create_gff)
        # The vfdb_results needs to be joined with locus tag info
        locus_tag_info = self.locus_df.select(["gene", "locus_tag", "contig", "start", "stop", "strand"])

        self.vfdb_results = self.vfdb_results.join(locus_tag_info, how="left", on="gene")

        # reorder: locus_tag first
        cols = self.vfdb_results.columns
        cols_reordered = ["locus_tag"] + [c for c in cols if c != "locus_tag"]
        self.vfdb_results = self.vfdb_results.select(cols_reordered)

        # keep only desired columns and save
        self.vfdb_results = self.vfdb_results.select([
            "contig", "locus_tag", "vfdb_hit", "vfdb_alnScore", "vfdb_seqIdentity", "start", "stop", "strand",
        ])
        self.vfdb_results = self.vfdb_results.rename({
            "locus_tag": "gene",
        })

        if self.mmseqs_flag is True:
            self.vfdb_results = self.vfdb_results.sort("start")
            self.vfdb_results.write_csv(
                os.path.join(self.out_dir, "top_hits_vfdb.tsv"), separator="\t"
            )

        ######################################
        ##### update card with locus tag #####
        ######################################
        self.card_results = self.card_results.join(locus_tag_info, how="left", on="gene")
        cols = self.card_results.columns
        cols_reordered = ["locus_tag"] + [c for c in cols if c != "locus_tag"]
        self.card_results = self.card_results.select(cols_reordered)

        self.card_results = self.card_results.select([
            "contig", "locus_tag", "CARD_hit", "CARD_alnScore", "CARD_seqIdentity", "start", "stop", "strand",
        ])
        self.card_results = self.card_results.rename({
            "locus_tag": "gene",
            "CARD_hit": "card_hit",
            "CARD_alnScore": "card_alnScore",
            "CARD_seqIdentity": "card_seqIdentity",
        })
        if self.mmseqs_flag is True:
            self.card_results = self.card_results.sort("start")
            self.card_results.write_csv(
                os.path.join(self.out_dir, "top_hits_card.tsv"), separator="\t"
            )

    def create_txt(self):
        """
        Creates the _cds_functions.tsv and _length_gc_cds_density.tsv outputs.

        Per-contig counts (CDS by PHROG category, tRNAs, CRISPRs, tmRNAs,
        VFDB, CARD) are computed in a single ``group_by`` pass over each
        source DataFrame, then assembled into the long-format output table.
        """
        # get contigs
        contigs = self.length_df["contig"].cast(pl.Utf8)
        contigs_list = contigs.to_list()
        self.merged_df = self.merged_df.with_columns(pl.col("contig").cast(pl.Utf8))

        col_list = [
            "contig", "Method", "Region", "start", "stop", "score", "strand", "frame", "attributes",
        ]

        if self.skip_extra_annotations is False:
            _empty_gff = pl.DataFrame({c: pl.Series([], dtype=pl.Utf8) for c in col_list})

            try:
                trna_df = pl.read_csv(
                    os.path.join(self.out_dir, "trnascan_out.gff"),
                    separator="\t",
                    has_header=False,
                    new_columns=col_list,
                    comment_prefix="#",
                    truncate_ragged_lines=True,
                    infer_schema=False,
                )
            except pl.exceptions.NoDataError:
                trna_df = _empty_gff
            trna_df = trna_df.filter(
                (pl.col("Region") == "tRNA") | (pl.col("Region") == "pseudogene")
            )

            try:
                crispr_df = pl.read_csv(
                    os.path.join(self.out_dir, self.prefix + "_minced.gff"),
                    separator="\t",
                    has_header=False,
                    new_columns=col_list,
                    comment_prefix="#",
                    infer_schema=False,
                )
            except pl.exceptions.NoDataError:
                crispr_df = _empty_gff

            try:
                tmrna_df = pl.read_csv(
                    os.path.join(self.out_dir, self.prefix + "_aragorn.gff"),
                    separator="\t",
                    has_header=False,
                    new_columns=col_list,
                    comment_prefix="#",
                    truncate_ragged_lines=True,
                    infer_schema=False,
                )
            except pl.exceptions.NoDataError:
                tmrna_df = _empty_gff

        # ─── Per-contig CDS stats: count + PHROG category breakdown + length sum
        # The 11 PHROG-category labels in display order.
        cds_categories = [
            ("connector",       "connector"),
            ("metabolism",      "DNA, RNA and nucleotide metabolism"),
            ("head",            "head and packaging"),
            ("integration",     "integration and excision"),
            ("lysis",           "lysis"),
            ("moron",           "moron, auxiliary metabolic gene and host takeover"),
            ("other",           "other"),
            ("tail",            "tail"),
            ("transcription",   "transcription regulation"),
            ("unknown",         "unknown function"),
        ]

        cds_only_df = self.merged_df.filter(pl.col("Region") == "CDS")
        if cds_only_df.height > 0:
            cds_only_df = cds_only_df.with_columns([
                # function name extracted from attributes (matches old splitn logic)
                pl.col("attributes")
                  .str.splitn(";function=", 2).struct.field("field_1")
                  .str.splitn(";product=", 2).struct.field("field_0")
                  .alias("_function"),
                # |start − stop| feeds the coding-density calculation
                (pl.col("start").cast(pl.Int64) - pl.col("stop").cast(pl.Int64))
                    .abs()
                    .alias("_cds_len"),
            ])
            cds_stats = cds_only_df.group_by("contig", maintain_order=True).agg(
                [pl.len().alias("CDS"), pl.col("_cds_len").sum().alias("_cds_len_sum")]
                + [
                    (pl.col("_function") == label).sum().alias(short)
                    for short, label in cds_categories
                ]
            )
        else:
            cds_stats = pl.DataFrame(schema={
                "contig": pl.Utf8, "CDS": pl.UInt32, "_cds_len_sum": pl.Int64,
                **{short: pl.UInt32 for short, _ in cds_categories},
            })

        cds_stats_dict = {
            row["contig"]: row for row in cds_stats.iter_rows(named=True)
        }

        # ─── Per-contig RNA / CRISPR / tmRNA counts via group_by
        def _counts_by_contig(src_df):
            if src_df.height == 0:
                return {}
            return dict(
                src_df.group_by("contig")
                      .agg(pl.len().alias("_n"))
                      .iter_rows()
            )

        if self.skip_extra_annotations is False:
            trna_counts   = _counts_by_contig(trna_df)
            crispr_counts = _counts_by_contig(crispr_df)
            tmrna_counts  = _counts_by_contig(tmrna_df)

        # ─── VFDB / CARD counts.  v1.9.1 used str.contains(contig) on the hit
        # row's contig column — preserved here in case the test suite ever
        # depends on substring matching (would need explicit handling if two
        # contig names were substrings of one another).
        def _counts_by_contig_contains(src_df, contig_names):
            if src_df.height == 0:
                return {c: 0 for c in contig_names}
            return {
                c: src_df.filter(pl.col("contig").str.contains(c)).height
                for c in contig_names
            }

        vfdb_counts = _counts_by_contig_contains(self.vfdb_results, contigs_list)
        card_counts = _counts_by_contig_contains(self.card_results, contigs_list)

        # ─── Total contig length for coding-density calculation
        length_lookup = dict(
            self.length_df.select(["contig", "length"]).iter_rows()
        )

        # ─── Assemble the long-format output table.  Build flat lists once,
        # then construct a single DataFrame (avoids 16+ tiny DataFrames/contig).
        descriptions, counts, contigs_out = [], [], []
        density_dict = {}

        for contig in contigs_list:
            stats = cds_stats_dict.get(contig)
            if stats is None:
                cds_count, cds_len_sum = 0, 0
                cat_counts = [0] * len(cds_categories)
            else:
                cds_count    = stats["CDS"]
                cds_len_sum  = stats["_cds_len_sum"]
                cat_counts   = [stats[short] for short, _ in cds_categories]

            # 11 PHROG-category rows
            descriptions.append("CDS"); counts.append(cds_count); contigs_out.append(contig)
            for (_short, label), n in zip(cds_categories, cat_counts):
                descriptions.append(label); counts.append(n); contigs_out.append(contig)

            if self.skip_extra_annotations is False:
                descriptions.append("tRNAs");   counts.append(trna_counts.get(contig, 0));   contigs_out.append(contig)
                descriptions.append("CRISPRs"); counts.append(crispr_counts.get(contig, 0)); contigs_out.append(contig)
                descriptions.append("tmRNAs");  counts.append(tmrna_counts.get(contig, 0));  contigs_out.append(contig)

            descriptions.append("VFDB_Virulence_Factors"); counts.append(vfdb_counts.get(contig, 0)); contigs_out.append(contig)
            descriptions.append("CARD_AMR_Genes");         counts.append(card_counts.get(contig, 0)); contigs_out.append(contig)

            # coding density
            contig_length_val = length_lookup.get(contig)
            density_dict[contig] = (
                round(cds_len_sum * 100 / contig_length_val, 2)
                if contig_length_val else 0
            )

        # update length_df with densities
        density_vals = [density_dict.get(c, None) for c in self.length_df["contig"].to_list()]
        self.length_df = self.length_df.with_columns(
            pl.Series("cds_coding_density", density_vals)
        )

        description_total_df = pl.DataFrame({
            "Description": descriptions,
            "Count": counts,
            "contig": contigs_out,
        })
        description_total_df.write_csv(
            os.path.join(self.out_dir, self.prefix + "_cds_functions.tsv"),
            separator="\t",
        )
        self.length_df.write_csv(
            os.path.join(self.out_dir, self.prefix + "_length_gc_cds_density.tsv"),
            separator="\t",
        )

    def update_fasta_headers(self):
        """
        Updates the fasta output headers to have a consistent locus tag & gene description
        """

        fasta_input_nts_tmp = self.gene_predictor + "_out_tmp.fasta"
        fasta_input_aas_tmp = self.gene_predictor + "_aas_tmp.fasta"
        fasta_output_nts_gd = self.gene_predictor + ".ffn"
        fasta_output_aas_gd = self.gene_predictor + ".faa"

        locus_tags = self.locus_df["locus_tag"].to_list()
        annots = self.locus_df["annot"].to_list()

        with open(os.path.join(self.out_dir, fasta_output_nts_gd), "w") as nt_fa:
            for i, dna_record in enumerate(SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_nts_tmp), "fasta"
            )):
                dna_record.id = str(locus_tags[i])
                dna_record.description = str(annots[i])
                SeqIO.write(dna_record, nt_fa, "fasta")

        with open(os.path.join(self.out_dir, fasta_output_aas_gd), "w") as aa_fa:
            for i, dna_record in enumerate(SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_aas_tmp), "fasta"
            )):
                dna_record.id = str(locus_tags[i])
                dna_record.description = str(annots[i])
                SeqIO.write(dna_record, aa_fa, "fasta")

    def update_final_output(self):
        """
        Updates the fasta output headers to have consistent locus tag & gene description
        """
        # rename gene with locus_tag
        locus_tag_list = self.locus_df["locus_tag"].to_list()
        self.merged_df = self.merged_df.with_columns(
            pl.Series("gene", locus_tag_list)
        )

        # rearrange start and stop for neg strand
        self.merged_df = self.merged_df.with_columns([
            pl.when(pl.col("strand") == "-").then(pl.col("stop")).otherwise(pl.col("start")).alias("start"),
            pl.when(pl.col("strand") == "-").then(pl.col("start")).otherwise(pl.col("stop")).alias("stop"),
        ])

        # move gene column to head
        cols = self.merged_df.columns
        cols_reordered = ["gene"] + [c for c in cols if c != "gene"]
        self.merged_df = self.merged_df.select(cols_reordered)

        # drop cols
        drop_cols = [c for c in ["frame", "attributes", "count", "locus_tag"] if c in self.merged_df.columns]
        self.merged_df = self.merged_df.drop(drop_cols)

        final_output_path = os.path.join(
            self.out_dir, self.prefix + "_cds_final_merged_output.tsv"
        )
        self.merged_df.write_csv(final_output_path, separator="\t")

    def extract_terl(self):
        """
        Extract large terminase subunit
        """

        fasta_input_nts_tmp = self.gene_predictor + "_out_tmp.fasta"
        fasta_input_aas_tmp = self.gene_predictor + "_aas_tmp.fasta"

        locus_tags = self.locus_df["locus_tag"].to_list()
        annots = self.locus_df["annot"].to_list()

        with open(os.path.join(self.out_dir, "terL.ffn"), "w") as aa_fa:
            j = 0
            for i, dna_record in enumerate(SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_nts_tmp), "fasta"
            )):
                dna_record.id = str(locus_tags[i])
                dna_record.description = str(annots[i])
                if annots[i] == "terminase large subunit":
                    SeqIO.write(dna_record, aa_fa, "fasta")
                    if j < 1:
                        logger.info("Terminase large subunit found.")
                    if j == 1:
                        logger.info(
                            "More than one CDS annotated as terminase large subunit found. \nSaving all."
                        )
                    j += 1

        with open(os.path.join(self.out_dir, "terL.faa"), "w") as aa_fa:
            for i, dna_record in enumerate(SeqIO.parse(
                os.path.join(self.out_dir, fasta_input_aas_tmp), "fasta"
            )):
                dna_record.id = str(locus_tags[i])
                dna_record.description = str(annots[i])
                if annots[i] == "terminase large subunit":
                    SeqIO.write(dna_record, aa_fa, "fasta")

    def inphared_top_hits(self):
        """
        Process mash output to get inphared top hits
        """

        mash_tsv = os.path.join(self.out_dir, "mash_out.tsv")
        col_list = [
            "contig", "Accession", "mash_distance", "mash_pval", "mash_matching_hashes",
        ]

        contigs = self.length_df["contig"].cast(pl.Utf8)
        contigs_list = contigs.to_list()

        try:
            mash_df = pl.read_csv(
                mash_tsv, separator="\t", has_header=False, new_columns=col_list
            )
            mash_df = mash_df.with_columns([
                pl.col("contig").cast(pl.Utf8),
                pl.col("mash_distance").cast(pl.Float64),
            ])
        except pl.exceptions.NoDataError:
            mash_df = pl.DataFrame({c: pl.Series([], dtype=pl.Utf8) for c in col_list})

        tophits = []

        for contig in contigs_list:
            # Sort by mash_distance only, preserving the original mash output row order for ties
            # (maintain_order=True is a stable sort).  Old pharokka used pandas sort_values
            # which is also stable, so ties were broken by the order in mash_out.tsv, which
            # reflects the INPHARED database row order (NOT alphabetical by Accession).
            # Differs from pharokka v1.9.1: previous polars code used a secondary Accession
            # sort which picked alphabetically-first accession for ties, e.g. "DQ222853"
            # instead of "NC_007458" (which appears earlier in INPHARED).
            # Example: DQ222853 (old, alphabetical) → NC_007458 (new, DB order for NC_007458.1 contig).
            hit_df = mash_df.filter(pl.col("contig") == contig).sort("mash_distance", maintain_order=True)
            hits = hit_df.height
            if hits > 0:
                top_row = hit_df[0]
                # mash_distance and mash_pval: always format as Python float strings.
                # Differs from pharokka v1.9.1 (pandas): pandas always read mash_pval
                # as float64, writing "0.0" for zero. Previous polars code depended on
                # dtype inference — writing "0" when polars inferred Int64 (all-zero
                # pval column) or "0.0" when it inferred Float64. New behaviour always
                # uses float() for a consistent representation.
                # Example for pval: "0" (old, dtype-dependent) → "0.0" (new, consistent).
                tophits.append([
                    top_row[0, "contig"],
                    top_row[0, "Accession"],
                    str(float(top_row[0, "mash_distance"])),
                    str(float(top_row[0, "mash_pval"])),
                    top_row[0, "mash_matching_hashes"],
                ])
            else:
                tophits.append([
                    contig,
                    "no_inphared_mash_hit",
                    "no_inphared_mash_hit",
                    "no_inphared_mash_hit",
                    "no_inphared_mash_hit",
                ])

        tophits_mash_df = pl.DataFrame(
            tophits,
            schema=["contig", "Accession", "mash_distance", "mash_pval", "mash_matching_hashes"],
            orient="row",
        )

        # read in the inphared tsv
        inphared_tsv_file = os.path.join(self.db_dir, "9Aug2025_data.tsv")
        cols = [
            "Accession", "Description", "Classification", "Genome_Length_(bp)", "Jumbophage",
            "molGC_(%)", "Molecule", "Modification_Date", "Number_CDS", "Positive_Strand_(%)",
            "Negative_Strand_(%)", "Coding_Capacity_(%)", "Low_Coding_Capacity_Warning",
            "tRNAs", "Host", "Lowest_Taxa", "Genus", "Sub-family", "Family", "Order",
            "Class", "Phylum", "Kingdom", "Realm", "Baltimore_Group", "Genbank_Division",
            "Isolation_Host_(beware_inconsistent_and_nonsense_values)",
        ]
        inphared_df = pl.read_csv(
            inphared_tsv_file,
            separator="\t",
            skip_rows=1,
            has_header=False,
            new_columns=cols,
            infer_schema=False,
        )

        # Normalise INPHARED columns to match old pandas behaviour:
        #   1. "NA" strings → null (pandas treated "NA" as NaN → wrote as empty, no quotes)
        #   2. "TRUE"/"FALSE" → "True"/"False" (pandas read as Python bool)
        #   3. Float cols with trailing zeros → stripped (pandas parsed as float64)
        inphared_df = inphared_df.with_columns([
            pl.when(pl.col(c) == "NA").then(None).otherwise(pl.col(c)).alias(c)
            for c in inphared_df.columns
        ])
        inphared_df = inphared_df.with_columns(
            pl.col("Jumbophage")
              .str.replace_all("^TRUE$", "True")
              .str.replace_all("^FALSE$", "False")
        )
        # Re-format float columns to strip trailing zeros (e.g. "44.890" → "44.89")
        float_inphared_cols = [
            "molGC_(%)", "Positive_Strand_(%)", "Negative_Strand_(%)", "Coding_Capacity_(%)",
        ]
        for _fc in float_inphared_cols:
            if _fc in inphared_df.columns:
                inphared_df = inphared_df.with_columns(
                    pl.when(pl.col(_fc).is_not_null())
                      .then(pl.col(_fc).cast(pl.Float64).cast(pl.Utf8))
                      .otherwise(pl.col(_fc))
                      .alias(_fc)
                )
        # NOTE (future cleanup): Genome_Length_(bp), Number_CDS, tRNAs were previously
        # cast to Float64 here to match an assumed pandas float64 behaviour.  However,
        # the 9Aug2025_data.tsv has zero NAs in those columns, so pandas reads them as
        # int64 and writes plain integers ("43411" not "43411.0").  With infer_schema=False
        # polars already preserves the raw integer strings, which is correct — no cast needed.

        combined_df = tophits_mash_df.join(inphared_df, on="Accession", how="left")
        combined_df.write_csv(
            os.path.join(self.out_dir, self.prefix + "_top_hits_mash_inphared.tsv"),
            separator="\t",
        )


#########################################
####### non class functions #############
########################################


def create_mmseqs_tophits(out_dir, reverse_mmseqs):
    """
    creates tophits_df dataframe from mmseqs2 phrog results
    """

    mmseqs_file = os.path.join(out_dir, "mmseqs_results.tsv")
    logger.info("Processing MMseqs2 outputs.")
    logger.info("Processing PHROGs output.")

    if reverse_mmseqs:
        col_list = [
            "gene", "mmseqs_phrog", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        ]
    else:
        col_list = [
            "mmseqs_phrog", "gene", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        ]

    try:
        mmseqs_df = pl.read_csv(
            mmseqs_file, separator="\t", has_header=False, new_columns=col_list, infer_schema=False
        )
    except pl.exceptions.NoDataError:
        mmseqs_df = pl.DataFrame({c: pl.Series([], dtype=pl.Utf8) for c in col_list})

    # Group by 'gene' and find the top hit for each group.
    # Sort by eVal (float) for correct ordering but preserve the original eVal
    # string so we can re-parse it with Python's float() later — polars uses a
    # different string→float64 parser than Python/numpy, causing tiny rounding
    # differences (e.g. "3.473e-22" → polars 3.4730000000000002e-22 vs Python
    # 3.4729999999999998e-22).
    tophits_df = (
        mmseqs_df
        .with_columns(pl.col("mmseqs_eVal").cast(pl.Float64).alias("_eVal_sort"))
        .sort("_eVal_sort", maintain_order=True)
        .unique(subset=["gene"], keep="first", maintain_order=True)
        .drop("_eVal_sort")
    )

    tophits_df = tophits_df.select([
        "mmseqs_phrog", "gene", "mmseqs_alnScore", "mmseqs_seqIdentity", "mmseqs_eVal",
    ])

    # Cast seqIdentity to float so trailing zeros are stripped ("0.660" → "0.66").
    # mmseqs_eVal: normalise to canonical float64 string via polars Ryu formatter.
    # Differs from pharokka v1.9.1 (pandas): mmseqs writes uppercase-E notation
    # (e.g. "3.473E-22"); pandas C-level float parser produced a different IEEE 754
    # bit pattern than polars for certain values, giving a longer repr.
    # Example: "3.4729999999999998e-22" (old, pandas parsed "3.473E-22") →
    #          "3.473e-22" (new, polars Ryu shortest round-trip). Also normalises
    # uppercase-E to lowercase-e.
    tophits_df = tophits_df.with_columns([
        pl.col("mmseqs_seqIdentity").cast(pl.Float64),
        pl.col("mmseqs_eVal").cast(pl.Float64).cast(pl.Utf8),
    ])

    tophits_df.write_csv(
        os.path.join(out_dir, "top_hits_mmseqs.tsv"), separator="\t"
    )
    return tophits_df


def remove_post_processing_files(out_dir, gene_predictor, meta, keep_raw_prodigal):
    """
    Cleans temporary files up
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

    remove_file(os.path.join(out_dir, "pharokka_tmp.gff"))
    remove_file(os.path.join(out_dir, "mash_out.tsv"))
    remove_file(os.path.join(out_dir, "input_mash_sketch.msh"))

    remove_file(os.path.join(out_dir, gene_predictor + "_aas_tmp.fasta"))

    if gene_predictor == "phanotate":
        remove_file(os.path.join(out_dir, "phanotate_out.txt"))
        remove_file(os.path.join(out_dir, gene_predictor + "_out_tmp.fasta"))
    elif gene_predictor == "genbank":
        remove_file(os.path.join(out_dir, "genbank.fasta"))
    else:  # prodigal or prodigal-gv
        remove_file(os.path.join(out_dir, gene_predictor + "_out.gff"))
        if keep_raw_prodigal:
            rename_file(
                os.path.join(out_dir, gene_predictor + "_out_tmp.fasta"),
                os.path.join(out_dir, gene_predictor + "_raw.ffn"),
            )
            rename_file(
                os.path.join(out_dir, gene_predictor + "_out_aas_tmp.fasta"),
                os.path.join(out_dir, gene_predictor + "_raw.faa"),
            )
        else:
            remove_file(os.path.join(out_dir, gene_predictor + "_out_tmp.fasta"))
            remove_file(os.path.join(out_dir, gene_predictor + "_out_aas_tmp.fasta"))

    if meta:
        remove_directory(os.path.join(out_dir, "input_split_tmp/"))


def get_crispr_count(out_dir, prefix):
    """
    Gets number of crisprs
    """
    crispr_file = os.path.join(out_dir, prefix + "_minced.gff")
    with open(crispr_file) as file:
        lines = file.readlines()
    crispr_count = 0
    for line in lines:
        if line[0] != "#":
            crispr_count += 1
    return crispr_count


def is_trna_empty(out_dir):
    """
    Determines if trna output file is empty
    """
    trna_empty = False
    if os.stat(os.path.join(out_dir, "trnascan_out.gff")).st_size == 0:
        trna_empty = True
    return trna_empty


def extract_anticodon_positions(out_dir):
    """Extract anticodon genomic positions from a tRNAscan-SE -f output file.

    Returns a dict mapping each tRNAscan tRNA ID (e.g. ``"MW460250_1.trna3"``)
    to its ``(anticodon_genomic_start, anticodon_genomic_end)`` tuple.

    The ``.sec`` file format alternates header lines like::

        MW460250_1.trna3 (115265-115194)\tLength: 72 bp
        Type: Met\tAnticodon: CAT at 33-35 (115233-115231)\tScore: 61.2

    Previously this function returned a positional list, which was matched
    against ``trna_df`` rows by index in ``create_gff``.  ``trna_df`` is
    sorted by genomic position (the GFF's order); the ``.sec`` file is
    sorted by tRNA ID — so for any contig with more than one tRNA, the
    anticodon ``(pos:X..Y)`` attribute in the output GFF/GBK was attached
    to the wrong tRNA.  Keying by the tRNAscan ID fixes that misalignment.
    """
    trna_ss_path = os.path.join(out_dir, "trnascan_out.sec")
    positions: dict = {}
    current_id = None
    with open(trna_ss_path) as f:
        for line in f:
            if line.startswith("Type:"):
                # Body line — attach the anticodon position to the ID we
                # parsed from the most recent header line.
                if current_id is None:
                    continue
                parts = line.split("Anticodon:")
                if len(parts) > 1:
                    after_anticodon = parts[1]
                    if "(" in after_anticodon and ")" in after_anticodon:
                        pos_str = after_anticodon.split("(")[1].split(")")[0]
                        start, end = pos_str.split("-")
                        positions[current_id] = (int(start), int(end))
                current_id = None
            else:
                # Header lines look like:
                #   "MW460250_1.trna3 (115265-115194)\tLength: 72 bp"
                # The tRNA ID is the first whitespace-delimited token and
                # contains ".trna" (distinguishes from blank / Seq: / Str: lines).
                token = line.split(maxsplit=1)[0] if line.strip() else ""
                if ".trna" in token:
                    current_id = token
    return positions


#### process pyhmmer hits
def process_pyhmmer_results(merged_df, pyhmmer_results_dict):
    """
    Processes pyhmmer
    """

    # split to get protein name to match with pyhmmer
    genes = merged_df["gene"].to_list()
    temp_prots = [g.split(" ", 1)[0] for g in genes]

    pyhmmer_phrog_list = ["No_PHROGs_HMM"] * len(genes)
    pyhmmer_bitscore_list = ["No_PHROGs_HMM"] * len(genes)
    pyhmmer_evalue_list = ["No_PHROGs_HMM"] * len(genes)

    for i, prot in enumerate(temp_prots):
        hit = pyhmmer_results_dict.get(prot)
        if hit is None:
            continue
        pyhmmer_phrog_list[i] = hit.phrog
        pyhmmer_bitscore_list[i] = str(round(hit.bitscore, 6))
        pyhmmer_evalue_list[i] = str(hit.evalue)

    merged_df = merged_df.with_columns([
        pl.Series("pyhmmer_phrog", pyhmmer_phrog_list),
        pl.Series("pyhmmer_bitscore", pyhmmer_bitscore_list),
        pl.Series("pyhmmer_evalue", pyhmmer_evalue_list),
    ])

    return merged_df


#### process pyhmmer hits
def process_custom_pyhmmer_results(merged_df, custom_pyhmmer_results_dict):
    """
    Processes pyhmmer results for custom db
    """

    genes = merged_df["gene"].to_list()
    temp_prots = [g.split(" ", 1)[0] for g in genes]

    custom_hmm_id_list = ["No_custom_HMM"] * len(genes)
    custom_hmm_bitscore_list = ["No_custom_HMM"] * len(genes)
    custom_hmm_evalue_list = ["No_custom_HMM"] * len(genes)

    for i, prot in enumerate(temp_prots):
        hit = custom_pyhmmer_results_dict.get(prot)
        if hit is None:
            continue
        custom_hmm_id_list[i] = hit.custom_hmm_id
        custom_hmm_bitscore_list[i] = str(round(hit.bitscore, 6))
        custom_hmm_evalue_list[i] = str(hit.evalue)

    merged_df = merged_df.with_columns([
        pl.Series("custom_hmm_id", custom_hmm_id_list),
        pl.Series("custom_hmm_bitscore", custom_hmm_bitscore_list),
        pl.Series("custom_hmm_evalue", custom_hmm_evalue_list),
    ])

    return merged_df


#### process vfdb files
def process_vfdb_results(out_dir, merged_df, proteins_flag=False, reverse_mmseqs2=False):
    """
    Processes VFDB results
    """
    vfdb_file = os.path.join(out_dir, "vfdb_results.tsv")
    logger.info("Processing VFDB output.")

    if reverse_mmseqs2:
        col_list = [
            "gene", "vfdb_hit", "vfdb_alnScore", "vfdb_seqIdentity", "vfdb_eVal",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        ]
    else:
        col_list = [
            "vfdb_hit", "gene", "vfdb_alnScore", "vfdb_seqIdentity", "vfdb_eVal",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        ]

    touch_file(vfdb_file)

    if os.path.getsize(vfdb_file) == 0:
        vfdb_df = pl.DataFrame({col: pl.Series([], dtype=pl.Utf8) for col in col_list})
    else:
        vfdb_df = pl.read_csv(
            vfdb_file, separator="\t", has_header=False, new_columns=col_list, infer_schema=False
        )

    # Issue #410 — strip problematic square brackets from certain VFDB hit
    # strings that break GenBank parsing.  All 8 substitutions are chained
    # on the same expression so polars processes them in a single pass over
    # the column rather than materialising 8 separate DataFrames.
    vfdb_df = vfdb_df.with_columns(
        pl.col("vfdb_hit")
          .str.replace_all(
              r"L-allo-isoleucine:holo-\[CmaA peptidyl-carrier protein\]",
              "L-allo-isoleucine:holo-CmaA peptidyl-carrier protein",
          )
          .str.replace_all(
              r"UDP-3-O-\[3-hydroxymyristoyl\]",
              "UDP-3-O-3-hydroxymyristoyl",
          )
          .str.replace_all(
              r"N-\[\(2S\)-2-amino-2-carboxyethyl\]",
              "N-(2S)-2-amino-2-carboxyethyl",
          )
          .str.replace_all(
              r"3-\(L-alanin-3-ylcarbamoyl\)-2-\[\(2- aminoethylcarbamoyl\)methyl\]",
              "3-(L-alanin-3-ylcarbamoyl)-2-(2- aminoethylcarbamoyl)methyl",
          )
          .str.replace_all(
              r"beta-ketoacyl-\[acyl-carrier-protein\]",
              "UDP-3-O-3-hydroxymyristoyl",
          )
          .str.replace_all(
              r"biotin--\[acetyl-CoA-carboxylase\] ligase",
              "biotin--acetyl-CoA-carboxylase ligase",
          )
          .str.replace_all(
              r"DP-3-O-\[3-hydroxymyristoyl\]",
              "DP-3-O-[3-hydroxymyristoyl]",
          )
          .str.replace_all(
              r"biotin--\[acetyl-CoA-carboxylase\]",
              "biotin--acetyl-CoA-carboxylase",
          )
    )

    # Sort by eVal; normalise eVal and seqIdentity to canonical float64 strings.
    tophits_df = (
        vfdb_df
        .with_columns(pl.col("vfdb_eVal").cast(pl.Float64).alias("_eVal_sort"))
        .sort("_eVal_sort", maintain_order=True)
        .unique(subset=["gene"], keep="first", maintain_order=True)
        .drop("_eVal_sort")
    )

    tophits_df = tophits_df.select([
        "vfdb_hit", "gene", "vfdb_alnScore", "vfdb_seqIdentity", "vfdb_eVal",
    ])

    # vfdb_eVal: normalise to canonical float64 string via polars Ryu formatter.
    # Differs from pharokka v1.9.1 (pandas): pandas C-level parser produced a
    # different IEEE 754 bit pattern for certain values, giving a longer repr.
    # Example: "2.807e-104" (old) vs "2.807e-104" (new, same here since it's
    # already the shortest round-trip). Exponent casing also normalised (E→e).
    #
    # vfdb_alnScore is the raw integer string from the mmseqs output (e.g. "321").
    # Differs from pharokka v1.9.1 (pandas): when null rows exist (CDS with no
    # VFDB hit), pandas promoted int64 to float64 due to NaN, writing "321.0".
    # New polars behaviour keeps the raw integer string.
    # Example: "321.0" (old, when ≥1 CDS had no VFDB hit) → "321" (new).
    tophits_df = tophits_df.with_columns([
        pl.col("vfdb_eVal").cast(pl.Float64).cast(pl.Utf8),
        pl.col("vfdb_seqIdentity").cast(pl.Float64).cast(pl.Utf8),
    ])

    tophits_df = tophits_df.with_columns(pl.col("gene").cast(pl.Utf8))

    if proteins_flag is True:
        tophits_df = tophits_df.with_columns(
            pl.col("gene").str.split(" ").list.first().alias("gene")
        )

    # merge top hits into the merged df
    merged_df = merged_df.join(tophits_df, on="gene", how="left")
    merged_df = merged_df.with_columns([
        pl.col("vfdb_hit").fill_null("None"),
        pl.col("vfdb_alnScore").fill_null("None"),
        pl.col("vfdb_seqIdentity").fill_null("None"),
        pl.col("vfdb_eVal").fill_null("None"),
    ])

    if tophits_df.height > 0:
        number_vfs = tophits_df.height
        logger.info(str(number_vfs) + " VFDB virulence factors identified.")

        merged_df = merged_df.with_columns(
            pl.col("vfdb_hit").str.splitn("[", 3).alias("_split")
        ).with_columns([
            pl.col("_split").struct.field("field_0").alias("genbank"),
            pl.col("_split").struct.field("field_1").alias("desc_tmp"),
            pl.col("_split").struct.field("field_2").alias("vfdb_species"),
        ]).drop("_split")

        merged_df = merged_df.with_columns(
            pl.col("vfdb_species").str.replace_all("]", "", literal=True).str.strip_chars().alias("vfdb_species")
        )

        # genbank has the info
        merged_df = merged_df.with_columns([
            pl.col("genbank").str.splitn(")", 2).struct.field("field_1").alias("vfdb_short_name"),
            pl.col("genbank").str.splitn(")", 3).struct.field("field_2").alias("vfdb_description"),
        ])
        merged_df = merged_df.with_columns(
            pl.col("vfdb_short_name").str.replace_all("(", "", literal=True)
            .str.splitn(")", 2).struct.field("field_0")
            .str.strip_chars()
            .alias("vfdb_short_name")
        )
        merged_df = merged_df.with_columns(
            pl.col("vfdb_description").str.strip_chars().alias("vfdb_description")
        )

        merged_df = merged_df.drop(["genbank", "desc_tmp"])
        merged_df = merged_df.with_columns([
            pl.col("vfdb_short_name").fill_null("None"),
            pl.col("vfdb_description").fill_null("None"),
            pl.col("vfdb_species").fill_null("None"),
        ])
    else:
        logger.info("0 VFDB virulence factors identified.")
        merged_df = merged_df.with_columns([
            pl.lit("None").alias("vfdb_short_name"),
            pl.lit("None").alias("vfdb_description"),
            pl.lit("None").alias("vfdb_species"),
        ])
    return (merged_df, tophits_df)


#### process CARD files
def process_card_results(out_dir, merged_df, db_dir, proteins_flag=False, reverse_mmseqs2=False):
    """
    Processes card results
    :param out_dir: output directory path
    :param merged_df: merged_df in process_results
    :proteins_flag bool, True if pharokka_proteins is run
    :return: merged_df merged_df updated with card results
    """
    card_file = os.path.join(out_dir, "CARD_results.tsv")
    logger.info("Processing CARD output.")
    if reverse_mmseqs2:
        col_list = [
            "gene", "CARD_hit", "CARD_alnScore", "CARD_seqIdentity", "CARD_eVal",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        ]
    else:
        col_list = [
            "CARD_hit", "gene", "CARD_alnScore", "CARD_seqIdentity", "CARD_eVal",
            "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen",
        ]
    touch_file(card_file)

    if os.path.getsize(card_file) == 0:
        card_df = pl.DataFrame({col: pl.Series([], dtype=pl.Utf8) for col in col_list})
    else:
        card_df = pl.read_csv(
            card_file, separator="\t", has_header=False, new_columns=col_list, infer_schema=False
        )

    # Sort by eVal; normalise eVal and seqIdentity to canonical float64 strings.
    tophits_df = (
        card_df
        .with_columns(pl.col("CARD_eVal").cast(pl.Float64).alias("_eVal_sort"))
        .sort("_eVal_sort", maintain_order=True)
        .unique(subset=["gene"], keep="first", maintain_order=True)
        .drop("_eVal_sort")
    )

    tophits_df = tophits_df.select([
        "CARD_hit", "gene", "CARD_alnScore", "CARD_seqIdentity", "CARD_eVal",
    ])

    # CARD_eVal: normalise to canonical float64 string via polars Ryu formatter.
    # Differs from pharokka v1.9.1 (pandas): pandas C-level parser produced a
    # different IEEE 754 bit pattern for certain values, giving a longer repr.
    # Example: same as mmseqs_eVal — "3.4729999999999998e-22" (old) →
    #          "3.473e-22" (new, polars Ryu shortest round-trip). Exponent E→e.
    #
    # CARD_alnScore is the raw integer string from the mmseqs output (e.g. "243").
    # Differs from pharokka v1.9.1 (pandas): when null rows exist (CDS with no
    # CARD hit), pandas promoted int64 to float64, writing "243.0".
    # New polars behaviour keeps the raw integer string.
    # Example: "243.0" (old, when ≥1 CDS had no CARD hit) → "243" (new).
    tophits_df = tophits_df.with_columns([
        pl.col("CARD_eVal").cast(pl.Float64).cast(pl.Utf8),
        pl.col("CARD_seqIdentity").cast(pl.Float64).cast(pl.Utf8),
    ])

    tophits_df = tophits_df.with_columns(pl.col("gene").cast(pl.Utf8))

    if proteins_flag is True:
        tophits_df = tophits_df.with_columns(
            pl.col("gene").str.split(" ").list.first().alias("gene")
        )

    merged_df = merged_df.join(tophits_df, on="gene", how="left")
    merged_df = merged_df.with_columns([
        pl.col("CARD_hit").fill_null("None"),
        pl.col("CARD_alnScore").fill_null("None"),
        pl.col("CARD_seqIdentity").fill_null("None"),
        pl.col("CARD_eVal").fill_null("None"),
    ])

    if tophits_df.height > 0:
        number_cards = tophits_df.height
        logger.info(str(number_cards) + " CARD AMR genes identified.")

        merged_df = merged_df.with_columns(
            pl.col("CARD_hit").str.splitn("[", 2).alias("_split")
        ).with_columns([
            pl.col("_split").struct.field("field_0").alias("genbank"),
            pl.col("_split").struct.field("field_1").alias("CARD_species"),
        ]).drop("_split")

        merged_df = merged_df.with_columns(
            pl.col("CARD_species").str.replace_all("]", "", literal=True).str.strip_chars().alias("CARD_species")
        )

        merged_df = merged_df.with_columns(
            pl.col("genbank").str.splitn("|", 4).alias("_split2")
        ).with_columns([
            pl.col("_split2").struct.field("field_0").alias("gb"),
            pl.col("_split2").struct.field("field_1").alias("genbank2"),
            pl.col("_split2").struct.field("field_2").alias("ARO_Accession"),
            pl.col("_split2").struct.field("field_3").alias("CARD_short_name"),
        ]).drop(["_split2", "genbank"]).rename({"genbank2": "genbank"})

        merged_df = merged_df.with_columns(
            pl.col("CARD_short_name").str.strip_chars().alias("CARD_short_name")
        )

        CARD_index_file = os.path.join(db_dir, "aro_index.tsv")
        card_cols = [
            "ARO_Accession", "CVTERM_ID", "Model_Sequence_ID", "Model_ID", "Model_Name",
            "ARO_Name", "Protein_Accession", "DNA_Accession", "AMR_Gene_Family",
            "Drug_Class", "Resistance_Mechanism", "CARD_Short_Name",
        ]
        card_index_df = pl.read_csv(
            CARD_index_file,
            separator="\t",
            skip_rows=1,
            has_header=False,
            new_columns=card_cols,
            infer_schema=False,
        )
        card_index_df = card_index_df.drop([
            "CVTERM_ID", "Model_Sequence_ID", "Model_ID", "Model_Name",
            "ARO_Name", "CARD_Short_Name",
        ])
        merged_df = merged_df.join(card_index_df, on="ARO_Accession", how="left")
        merged_df = merged_df.drop(["gb", "genbank"])

        merged_df = merged_df.with_columns([
            pl.col("CARD_species").fill_null("None"),
            pl.col("ARO_Accession").fill_null("None"),
            pl.col("CARD_short_name").fill_null("None"),
            pl.col("Protein_Accession").fill_null("None"),
            pl.col("DNA_Accession").fill_null("None"),
            pl.col("AMR_Gene_Family").fill_null("None"),
            pl.col("Drug_Class").fill_null("None"),
            pl.col("Resistance_Mechanism").fill_null("None"),
        ])
    else:
        logger.info("0 CARD AMR genes identified.")
        merged_df = merged_df.with_columns([
            pl.lit("None").alias("CARD_species"),
            pl.lit("None").alias("ARO_Accession"),
            pl.lit("None").alias("CARD_short_name"),
            pl.lit("None").alias("Protein_Accession"),
            pl.lit("None").alias("DNA_Accession"),
            pl.lit("None").alias("AMR_Gene_Family"),
            pl.lit("None").alias("Drug_Class"),
            pl.lit("None").alias("Resistance_Mechanism"),
        ])

    return (merged_df, tophits_df)


def is_file_empty(file):
    """
    Determines if file is empty
    """
    empty = False
    if os.stat(file).st_size == 0:
        empty = True
    return empty


def check_and_create_directory(directory):
    """Checks if directory exists, creates it if not"""
    if not os.path.isdir(directory):
        os.mkdir(directory)


def parse_attributes_column(df):
    """Parse a GFF-style ``attributes`` column into separate columns.

    Attributes look like ``ID=foo;phrog=123;product=bar``.  Each key=value
    pair becomes its own column.  Keys missing from a row produce nulls.

    Column ordering matches the order in which each key first appears across
    the rows (preserving the original loop-based behaviour).

    Vectorised implementation: splits the column on ``;`` and ``=`` using
    polars list/struct kernels, then pivots into wide form.  No Python loop
    over rows — runs in a single polars query for the whole column.
    """
    if df.height == 0 or "attributes" not in df.columns:
        return df

    # Tokenise: attributes → list[struct{key, value}], drop empty / malformed entries.
    # _pair_pos preserves the position of each pair within its row so that the
    # final column ordering matches the order in which the keys appear in the
    # original attributes string (matches the old loop-based dict iteration).
    pairs = (
        df.select(
            pl.int_range(0, df.height).alias("_row"),
            pl.col("attributes")
              .str.split(";")
              .list.eval(
                  pl.element().str.splitn("=", 2).struct.rename_fields(["key", "value"])
              )
              .alias("_pairs"),
        )
        .with_columns(
            pl.int_ranges(0, pl.col("_pairs").list.len()).alias("_pair_pos")
        )
        .explode(["_pairs", "_pair_pos"])
        .unnest("_pairs")
        .filter(pl.col("key").is_not_null() & (pl.col("key") != "") & pl.col("value").is_not_null())
    )

    if pairs.height == 0:
        return df

    # Preserve first-seen order of keys: within a single row, dict insertion
    # order = pair position; across rows, the first row each key appears in
    # wins.  Sort by (row, pair_pos) of first appearance.
    key_order = (
        pairs.group_by("key", maintain_order=True)
             .agg([
                 pl.col("_row").min().alias("_first_row"),
                 pl.col("_pair_pos").first().alias("_first_pos"),
             ])
             .sort(["_first_row", "_first_pos"])
             ["key"]
             .to_list()
    )

    # Pivot to wide form: one column per key, row-indexed by _row.
    wide = (
        pairs.pivot(
            on="key",
            index="_row",
            values="value",
            aggregate_function="first",
        )
        .sort("_row")
        .drop("_row")
    )

    # Reorder columns to match first-appearance order.
    wide = wide.select(key_order)

    cols_to_drop = [c for c in df.columns if c in wide.columns]
    return pl.concat([df.drop(cols_to_drop), wide], how="horizontal")


def get_codon_from_anticodon(anticodon):
    """
    Given a DNA sequence (string), return:
    - its reverse complement in RNA format

    Example:
        >>> get_codon_from_anticodon("ATG")
        'rna_from_revcomp': 'CAU'
    """
    if not isinstance(anticodon, str):
        raise TypeError("Input must be a string.")
    if not anticodon.strip():
        raise ValueError("Input DNA sequence is empty.")

    seq = Seq(anticodon.upper())

    revcomp = seq.reverse_complement()
    rna_from_revcomp = revcomp.transcribe()

    return str(rna_from_revcomp)
