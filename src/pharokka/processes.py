import multiprocessing.pool
import os
import subprocess as sp
from datetime import datetime

import polars as pl
import pyrodigal
import pyrodigal_gv
import pyrodigal_rv
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger

from .external_tools import ExternalTool
from .util import parse_attributes, remove_directory


def run_pyrodigal_gv(filepath_in, out_dir, threads):
    """
    Gets CDS using pyrodigal_gv
    """
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
                        genes.write_translations(aa_fasta, sequence_id=record_id)


def run_pyrodigal_rv(filepath_in, out_dir, threads):
    """
    Gets CDS using pyrodigal_rv
    """
    orf_finder = pyrodigal_rv.ViralGeneFinder(meta=True)

    def _find_genes(record):
        genes = orf_finder.find_genes(str(record.seq))
        return (record.id, genes)

    with multiprocessing.pool.ThreadPool(threads) as pool:
        with open(os.path.join(out_dir, "pyrodigal-rv_out.gff"), "w") as gff:
            with open(os.path.join(out_dir, "pyrodigal-rv_out_tmp.fasta"), "w") as dst:
                with open(
                    os.path.join(out_dir, "pyrodigal-rv_out_aas_tmp.fasta"), "w"
                ) as aa_fasta:
                    records = SeqIO.parse(filepath_in, "fasta")
                    for record_id, genes in pool.imap(_find_genes, records):
                        genes.write_gff(
                            gff, sequence_id=record_id, include_translation_table=True
                        )
                        genes.write_genes(dst, sequence_id=record_id)
                        genes.write_translations(aa_fasta, sequence_id=record_id)


##### phanotate meta mode ########


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size."""
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []


def split_input_fasta(filepath_in, out_dir):
    """Splits the input fasta into separate single fasta files for multithreading with phanotate."""
    input_tmp_dir = os.path.join(out_dir, "input_split_tmp")
    num_fastas = 0
    with open(filepath_in) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            filename = "input_subprocess%d.fasta" % (i + 1)
            with open(os.path.join(input_tmp_dir, filename), "w") as out:
                SeqIO.write([record], out, "fasta")
            num_fastas += 1
    return num_fastas


def run_phanotate_fasta_meta(filepath_in, out_dir, threads, num_fastas):
    """Runs phanotate to output fastas."""
    phanotate_tmp_dir = os.path.join(out_dir, "input_split_tmp")
    commands = []

    for i in range(1, num_fastas + 1):
        in_path = os.path.join(phanotate_tmp_dir, "input_subprocess" + str(i) + ".fasta")
        out_path = os.path.join(phanotate_tmp_dir, "phanotate_out_tmp" + str(i) + ".fasta")
        cmd = ["phanotate.py", in_path, "-o", out_path, "-f", "fasta"]
        commands.append(cmd)

    n = int(threads)
    for j in range(max(int(len(commands) / n) + 1, 1)):
        procs = [
            sp.Popen(cmd)
            for cmd in commands[j * n : min((j + 1) * n, len(commands))]
        ]
        for p in procs:
            p.wait()


def run_phanotate_txt_meta(filepath_in, out_dir, threads, num_fastas):
    """Runs phanotate to output text file."""
    phanotate_tmp_dir = os.path.join(out_dir, "input_split_tmp")
    commands = []

    for i in range(1, num_fastas + 1):
        in_path = os.path.join(phanotate_tmp_dir, "input_subprocess" + str(i) + ".fasta")
        out_path = os.path.join(phanotate_tmp_dir, "phanotate_out_tmp" + str(i) + ".txt")
        cmd = ["phanotate.py", in_path, "-o", out_path, "-f", "tabular"]
        commands.append(cmd)

    n = int(threads)
    for j in range(max(int(len(commands) / n) + 1, 1)):
        procs = [
            sp.Popen(cmd)
            for cmd in commands[j * n : min((j + 1) * n, len(commands))]
        ]
        for p in procs:
            p.wait()


def concat_phanotate_meta(out_dir, num_fastas):
    """Concatenates phanotate output for downstream analysis."""
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


def run_trnascan_meta(filepath_in, out_dir, threads, num_fastas, trna_scan_model):
    """Runs trnascan to output gffs one contig per thread."""
    input_tmp_dir = os.path.join(out_dir, "input_split_tmp")
    commands = []

    if trna_scan_model == "general":
        model = "G"
    else:
        model = "B"

    for i in range(1, num_fastas + 1):
        in_path = os.path.join(input_tmp_dir, "input_subprocess" + str(i) + ".fasta")
        filepath_out = os.path.join(input_tmp_dir, "trnascan_tmp" + str(i) + ".gff")
        sec_out = os.path.join(input_tmp_dir, "trnascan_tmp" + str(i) + ".sec")
        cmd = ["tRNAscan-SE", in_path, "--thread", "1", f"-{model}", "-D", "-Q", "-j", filepath_out, "-f", sec_out]
        commands.append(cmd)

    n = int(threads)
    for j in range(max(int(len(commands) / n) + 1, 1)):
        procs = [
            sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.DEVNULL)
            for cmd in commands[j * n : min((j + 1) * n, len(commands))]
        ]
        for p in procs:
            p.wait()


def concat_trnascan_meta(out_dir, num_fastas):
    """Concatenates trnascan output for downstream analysis."""
    input_tmp_dir = os.path.join(out_dir, "input_split_tmp")

    gffs = []
    for i in range(1, int(num_fastas) + 1):
        out_gff = "trnascan_tmp" + str(i) + ".gff"
        gffs.append(os.path.join(input_tmp_dir, out_gff))

    with open(os.path.join(out_dir, "trnascan_out.gff"), "w") as outfile:
        for fname in gffs:
            with open(fname) as infile:
                outfile.write(infile.read())

    sec_strucs = []
    for i in range(1, int(num_fastas) + 1):
        sec_struc_file = "trnascan_tmp" + str(i) + ".sec"
        sec_strucs.append(os.path.join(input_tmp_dir, sec_struc_file))

    with open(os.path.join(out_dir, "trnascan_out.sec"), "w") as outfile:
        for fname in sec_strucs:
            with open(fname) as infile:
                outfile.write(infile.read())


##### single contig mode ######


def run_phanotate(filepath_in, out_dir, logdir):
    """Runs phanotate."""
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
    except Exception:
        logger.error("Error with Phanotate\n")


def run_pyrodigal(filepath_in, out_dir, meta, coding_table, threads):
    """Gets CDS using pyrodigal."""
    prodigal_metamode = False
    if meta == True:
        prodigal_metamode = True
        logger.info("Prodigal Meta Mode Enabled")

    seqs = []
    total_length = 0
    with open(filepath_in, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs.append(bytes(record.seq))
            total_length += len(record.seq)

    if total_length < 100001:
        orf_finder = pyrodigal.GeneFinder(meta=True)
        prodigal_metamode = True
    else:
        orf_finder = pyrodigal.GeneFinder(meta=prodigal_metamode)
        if not prodigal_metamode:
            orf_finder.train(*seqs, translation_table=int(coding_table))

    def _find_genes(record):
        genes = orf_finder.find_genes(str(record.seq))
        return (record.id, genes)

    if prodigal_metamode:
        transl_table = None
    else:
        transl_table = "11"

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
                        genes.write_translations(aa_fasta, sequence_id=record_id)

                        if transl_table is None:
                            for gene in genes:
                                transl_table = str(gene.translation_table)
                                break
    return str(transl_table)


def tidy_phanotate_output(out_dir):
    """Tidies phanotate output — returns a polars DataFrame."""
    phan_file = os.path.join(out_dir, "phanotate_out.txt")
    col_list = ["start", "stop", "strand", "contig", "score"]

    phan_df = pl.read_csv(
        phan_file,
        separator="\t",
        has_header=False,
        new_columns=col_list,
        comment_prefix="#",
        schema_overrides={
            "start": pl.Utf8,
            "stop": pl.Utf8,
            "strand": pl.Utf8,
            "contig": pl.Utf8,
            "score": pl.Utf8,
        },
        infer_schema=False,
    )

    # Cast to proper types
    phan_df = phan_df.with_columns([
        pl.col("start").cast(pl.Int64),
        pl.col("stop").cast(pl.Int64),
    ])

    # Parse PHANOTATE score strings (E-notation, e.g. "-1.150410E+28") to
    # Float64. Polars uses a consistent IEEE 754 double representation for both
    # parsing and CSV output — no pandas/C-strtod quirks needed.
    # Old pharokka used pd.read_csv(dtype={"score": float}) then to_csv, which
    # applies C strtod (not Python float()) to parse phanotate's E-notation strings.
    # C strtod can give different IEEE 754 bit patterns for strings with trailing
    # zeros (e.g. "-1.150410E+28" → same bits as Python float → to_csv gives
    # "-1.15041e+28") versus those without (e.g. "-6.018731E+35" → different bits
    # → to_csv gives "-6.0187309999999995e+35").
    phan_df = phan_df.with_columns(pl.col("score").cast(pl.Float64))

    # Build gene column using row index
    phan_df = phan_df.with_row_index("_idx")
    phan_df = phan_df.with_columns(
        (
            pl.col("contig").cast(pl.Utf8)
            + pl.col("_idx").cast(pl.Utf8)
            + pl.lit(" ")
            + pl.col("start").cast(pl.Utf8)
            + pl.lit("_")
            + pl.col("stop").cast(pl.Utf8)
        ).alias("gene")
    ).drop("_idx")

    phan_df.write_csv(
        os.path.join(out_dir, "cleaned_phanotate.tsv"), separator="\t"
    )
    return phan_df


def tidy_prodigal_output(out_dir, gene_predictor):
    """Tidies prodigal output — returns a polars DataFrame."""
    prefix = gene_predictor
    prod_file = os.path.join(out_dir, f"{prefix}_out.gff")
    col_list = [
        "contig", "prod", "orf", "start", "stop",
        "score", "strand", "frame", "description",
    ]

    prod_df = pl.read_csv(
        prod_file,
        separator="\t",
        has_header=False,
        new_columns=col_list,
        comment_prefix="#",
        schema_overrides={
            "contig": pl.Utf8,
            "prod": pl.Utf8,
            "orf": pl.Utf8,
            "start": pl.Int64,
            "stop": pl.Int64,
            "score": pl.Float64,
            "strand": pl.Utf8,
            "frame": pl.Utf8,
            "description": pl.Utf8,
        },
        infer_schema=False,
    )

    # Drop nulls (meta mode can introduce them)
    prod_df = prod_df.drop_nulls()

    # Keep relevant columns
    prod_filt_df = prod_df.select(["start", "stop", "strand", "contig", "score", "description"])

    # Parse the description column to get attributes
    parsed_attrs = [parse_attributes(row) for row in prod_filt_df["description"].to_list()]

    # Build a polars DataFrame from the parsed attributes
    attr_keys = set()
    for d in parsed_attrs:
        attr_keys.update(d.keys())

    attr_data = {k: [d.get(k) for d in parsed_attrs] for k in attr_keys}
    # Remove "score" to avoid duplication
    attr_data.pop("score", None)
    attr_df = pl.DataFrame(attr_data)

    prod_filt_df = pl.concat([prod_filt_df, attr_df], how="horizontal")

    # Keep only needed columns (partial may not always exist)
    keep_cols = ["start", "stop", "strand", "contig", "score", "partial"]
    available_cols = [c for c in keep_cols if c in prod_filt_df.columns]
    prod_filt_df = prod_filt_df.select(available_cols)
    if "partial" not in prod_filt_df.columns:
        prod_filt_df = prod_filt_df.with_columns(pl.lit(None).cast(pl.Utf8).alias("partial"))

    prod_filt_df = prod_filt_df.with_columns([
        pl.col("start").cast(pl.Int64),
        pl.col("stop").cast(pl.Int64),
    ])

    # Swap start/stop for negative-strand genes (to match phanotate_out convention)
    prod_filt_df = prod_filt_df.with_columns([
        pl.when(pl.col("strand") == "-")
        .then(pl.col("stop"))
        .otherwise(pl.col("start"))
        .alias("start"),
        pl.when(pl.col("strand") == "-")
        .then(pl.col("start"))
        .otherwise(pl.col("stop"))
        .alias("stop"),
    ])

    prod_filt_df = prod_filt_df.with_row_index("_idx")
    prod_filt_df = prod_filt_df.with_columns(
        (
            pl.col("contig").cast(pl.Utf8)
            + pl.col("_idx").cast(pl.Utf8)
            + pl.lit(" ")
            + pl.col("start").cast(pl.Utf8)
            + pl.lit("_")
            + pl.col("stop").cast(pl.Utf8)
        ).alias("gene")
    ).drop("_idx")

    prod_filt_df.write_csv(
        os.path.join(out_dir, f"cleaned_{prefix}.tsv"), separator="\t"
    )
    return prod_filt_df


def tidy_genbank_output(out_dir, genbank_file, coding_table):
    """Tidies --genbank input output — returns a polars DataFrame."""
    starts = []
    stops = []
    strands = []
    contigs = []

    fasta_nucl_tmp = "genbank_out_tmp.fasta"
    fasta_output_aas_tmp = "genbank_aas_tmp.fasta"

    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                strand = feature.location.strand
                if strand == 1:
                    strand = "+"
                    start = feature.location.start + 1
                    stop = feature.location.end
                else:
                    strand = "-"
                    start = feature.location.end
                    stop = feature.location.start + 1
                contig = record.id
                starts.append(start)
                stops.append(stop)
                strands.append(strand)
                contigs.append(contig)

    gen_df = pl.DataFrame({
        "start": starts,
        "stop": stops,
        "strand": strands,
        "contig": contigs,
    })
    gen_df = gen_df.with_columns(pl.lit("No_score").alias("score"))

    gen_df = gen_df.with_row_index("_idx")
    gen_df = gen_df.with_columns(
        (
            pl.col("contig").cast(pl.Utf8)
            + pl.col("_idx").cast(pl.Utf8)
            + pl.lit(" ")
            + pl.col("start").cast(pl.Utf8)
            + pl.lit("_")
            + pl.col("stop").cast(pl.Utf8)
        ).alias("gene")
    ).drop("_idx")

    gen_df.write_csv(os.path.join(out_dir, "cleaned_genbank.tsv"), separator="\t")

    # Get list values for random-access in the loop below
    contig_list = gen_df["contig"].to_list()
    start_list = gen_df["start"].to_list()
    stop_list = gen_df["stop"].to_list()

    with open(os.path.join(out_dir, fasta_nucl_tmp), "w") as nt_fa:
        with open(os.path.join(out_dir, fasta_output_aas_tmp), "w") as aa_fa:
            i = 0
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        header = str(contig_list[i]) + str(i)
                        description = str(start_list[i]) + "_" + str(stop_list[i])

                        cds_nt_seq = feature.location.extract(record.seq)
                        cds_amino_acid_seq = cds_nt_seq.translate(
                            to_stop=True, table=coding_table
                        )

                        nt_record = SeqRecord(
                            cds_nt_seq,
                            id=header,
                            description=description,
                        )
                        aa_record = SeqRecord(
                            cds_amino_acid_seq,
                            id=header,
                            description=description,
                        )

                        SeqIO.write(aa_record, aa_fa, "fasta")
                        SeqIO.write(nt_record, nt_fa, "fasta")
                        i += 1

    return gen_df


def translate_fastas(out_dir, gene_predictor, coding_table, genbank_file):
    """Translates input CDSs to amino acids."""
    if gene_predictor == "phanotate":
        clean_df = tidy_phanotate_output(out_dir)
    elif gene_predictor in ("prodigal", "prodigal-gv", "pyrodigal-rv"):
        clean_df = tidy_prodigal_output(out_dir, gene_predictor)
    elif gene_predictor == "genbank":
        clean_df = tidy_genbank_output(out_dir, genbank_file, coding_table)

    fasta_output_aas_tmp = gene_predictor + "_aas_tmp.fasta"

    contig_list = clean_df["contig"].to_list()
    start_list = clean_df["start"].to_list()
    stop_list = clean_df["stop"].to_list()

    if gene_predictor == "phanotate":
        fasta_input_tmp = gene_predictor + "_out_tmp.fasta"
        with open(os.path.join(out_dir, fasta_output_aas_tmp), "w") as aa_fa:
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(out_dir, fasta_input_tmp), "fasta"
            ):
                dna_header = str(contig_list[i]) + str(i)
                dna_description = str(start_list[i]) + "_" + str(stop_list[i])
                aa_record = SeqRecord(
                    dna_record.seq.translate(to_stop=True, table=coding_table),
                    id=dna_header,
                    description=dna_description,
                )
                SeqIO.write(aa_record, aa_fa, "fasta")
                i += 1
    elif gene_predictor in ("prodigal-gv", "prodigal", "pyrodigal-rv"):
        fasta_input_tmp = gene_predictor + "_out_aas_tmp.fasta"
        with open(os.path.join(out_dir, fasta_output_aas_tmp), "w") as aa_fa:
            i = 0
            for dna_record in SeqIO.parse(
                os.path.join(out_dir, fasta_input_tmp), "fasta"
            ):
                dna_header = str(contig_list[i]) + str(i)
                dna_description = str(start_list[i]) + "_" + str(stop_list[i])
                aa_record = SeqRecord(
                    dna_record.seq,
                    id=dna_header,
                    description=dna_description,
                )
                SeqIO.write(aa_record, aa_fa, "fasta")
                i += 1
    # for genbank do nothing


def run_trna_scan(filepath_in, threads, out_dir, logdir, trna_scan_model):
    """Runs trna scan."""
    out_gff = os.path.join(out_dir, "trnascan_out.gff")
    out_sec = os.path.join(out_dir, "trnascan_out.sec")

    if trna_scan_model == "general":
        model = "G"
    else:
        model = "B"

    trna = ExternalTool(
        tool="tRNAscan-SE",
        input=f"{filepath_in}",
        output=f"{out_gff}",
        params=f"--thread {threads} -{model} -f {out_sec} -Q -j",
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(trna)
    except Exception:
        logger.error("Error with tRNAscan-SE")
        return 0


def run_mmseqs(db_dir, out_dir, threads, logdir, gene_predictor, evalue, reverse_mmseqs2, sensitivity, db_name):
    """Runs mmseqs2."""
    logger.info(f"Running MMseqs2 on {db_name} Database.")

    amino_acid_fasta = gene_predictor + "_aas_tmp.fasta"

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

    if os.path.isdir(target_db_dir) == False:
        os.mkdir(target_db_dir)

    mmseqs_createdb = ExternalTool(
        tool="mmseqs createdb",
        input=f"",
        output=f"{target_seqs}",
        params=f"{input_aa_fasta}",
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


def convert_gff_to_gbk(filepath_in, input_dir, out_dir, prefix, prot_seq_df):
    """Converts the gff to genbank.

    prot_seq_df is a polars DataFrame with columns 'contig' and 'sequence'.
    """
    gff_file = os.path.join(input_dir, f"{prefix}.gff")
    gbk_file = os.path.join(out_dir, f"{prefix}.gbk")

    prot_seq_df = prot_seq_df.with_columns(pl.col("contig").cast(pl.Utf8))

    with open(gbk_file, "wt") as gbk_handler:
        fasta_handler = SeqIO.to_dict(SeqIO.parse(filepath_in, "fasta"))

        for record in GFF.parse(gff_file, fasta_handler):
            record.id = str(record.id)
            sequence_length = len(record.seq)
            subset_df = prot_seq_df.filter(pl.col("contig") == record.id)
            subset_seqs = subset_df["sequence"].to_list()
            i = 0

            record.annotations["molecule_type"] = "DNA"
            record.annotations["date"] = datetime.today()
            record.annotations["topology"] = "linear"
            record.annotations["data_file_division"] = "PHG"

            for feature in record.features:
                if "source" in feature.qualifiers:
                    feature.qualifiers["inference"] = feature.qualifiers["source"]
                    del feature.qualifiers["source"]

                if "anticodon" in feature.qualifiers:
                    if isinstance(feature.qualifiers["anticodon"], list):
                        feature.qualifiers["anticodon"] = [
                            ",".join(feature.qualifiers["anticodon"])
                        ]

                if feature.type == "CDS":
                    if "phase" in feature.qualifiers:
                        del feature.qualifiers["phase"]

                    if "partial" in feature.qualifiers:
                        if (
                            feature.qualifiers["partial"] == ["10"]
                            and feature.location.strand == 1
                            and feature.location.start != 0
                        ):
                            feature.qualifiers["codon_start"] = (
                                feature.location.start + 1
                            )
                        elif (
                            feature.qualifiers["partial"] == ["01"]
                            and feature.location.strand == -1
                        ):
                            feature.qualifiers["codon_start"] = (
                                sequence_length - feature.location.end + 1
                            )

                    feature.qualifiers.update({"translation": subset_seqs[i]})
                    i += 1

            SeqIO.write(record, gbk_handler, "genbank")


def run_minced(filepath_in, out_dir, prefix, minced_args, logdir):
    """Runs MinCED."""
    logger.info("Running MinCED.")

    output_spacers = os.path.join(out_dir, prefix + "_minced_spacers.txt")
    output_gff = os.path.join(out_dir, prefix + "_minced.gff")

    if minced_args != "":
        minced_args = f"-{minced_args}"

    minced_fast = ExternalTool(
        tool="minced",
        input=f"",
        output=f" {output_spacers} {output_gff}",
        params=f" {minced_args} {filepath_in}",
        logdir=logdir,
        outfile="",
    )

    try:
        ExternalTool.run_tool(minced_fast)
    except Exception:
        logger.error("Error with MinCED\n")


def run_aragorn(filepath_in, out_dir, prefix, logdir):
    """Runs aragorn."""
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
    except Exception:
        logger.error("Error with Aragorn\n")


def reorient_terminase(filepath_in, out_dir, prefix, terminase_strand, terminase_start):
    """Re-orients phage to begin with large terminase subunit."""
    logger.info(
        "Input checked. \nReorienting input genome to begin with terminase large subunit."
    )

    record = SeqIO.read(filepath_in, "fasta")
    length = len(record.seq)

    if int(terminase_start) > length or int(terminase_start) < 1:
        logger.error(
            "Error: terminase large subunit start coordinate specified is not within the provided genome length. Please check your input. \n"
        )

    if terminase_strand == "pos":
        start = record.seq[(int(terminase_start) - 1) : length]
        end = record.seq[0 : int(terminase_start) - 1]
        total = start + end

    if terminase_strand == "neg":
        record.seq = record.seq.reverse_complement()
        start = record.seq[(length - int(terminase_start)) : length]
        end = record.seq[0 : (length - int(terminase_start))]
        total = start + end

    record.seq = total
    out_fasta = os.path.join(out_dir, prefix + "_genome_terminase_reoriented.fasta")
    SeqIO.write(record, out_fasta, "fasta")


def run_mash_sketch(filepath_in, out_dir, logdir):
    """Runs mash sketch."""
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
    except Exception:
        logger.error("Error with mash sketch\n")


def run_mash_dist(out_dir, db_dir, mash_distance, logdir):
    """Runs mash dist."""
    mash_sketch = os.path.join(out_dir, "input_mash_sketch.msh")
    phrog_sketch = os.path.join(db_dir, "9Aug2025_genomes.fa.msh")
    mash_tsv = os.path.join(out_dir, "mash_out.tsv")

    mash_dist = ExternalTool(
        tool="mash",
        input="",
        output="",
        params=f" dist  {mash_sketch} {phrog_sketch} -d {mash_distance} -i ",
        logdir=logdir,
        outfile=mash_tsv,
    )

    try:
        ExternalTool.run_tool(mash_dist, to_stdout=True)
    except Exception:
        logger.error("Error with mash dist\n")


def run_dnaapler(filepath_in, contig_count, out_dir, threads, logdir):
    """Runs dnaapler."""
    logger.info(
        "Running Dnaapler to rerorient all contigs to begin with the terminase large subunit."
    )

    dnaapler_output_dir = os.path.join(out_dir, "dnaapler")

    if contig_count == 1:
        dnaapler = ExternalTool(
            tool="dnaapler phage",
            input=f"-i {filepath_in}",
            output=f"-o {dnaapler_output_dir} -t {threads}",
            params=f"",
            logdir=logdir,
            outfile="",
        )
    else:
        dnaapler = ExternalTool(
            tool="dnaapler all",
            input=f"-i {filepath_in}",
            output=f"-o {dnaapler_output_dir} -t {threads}",
            params=f"",
            logdir=logdir,
            outfile="",
        )

    dnaapler_success = True
    try:
        ExternalTool.run_tool(dnaapler)
    except Exception:
        logger.warning("Dnaapler failed to find the large terminase subunit.")
        logger.warning(
            "For reorientation, please run pharokka again without --dnaapler and use --terminase instead."
        )
        logger.warning("Pharokka will now continue without reorientation.")
        dnaapler_success = False

    return dnaapler_success
