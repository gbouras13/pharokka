"""Parsing of Infernal cmscan output against Rfam.

The tabular output produced by ``cmscan --fmt 2`` is whitespace-delimited with
a fixed column order and a free-text description in the final field.  The
columns are (1-indexed, as documented in the Infernal user guide):

     1 idx              11 seq to          21 anyidx
     2 target name      12 strand          22 afrct1
     3 accession        13 trunc           23 afrct2
     4 query name       14 pass            24 winidx
     5 accession        15 gc              25 wfrct1
     6 clan name        16 bias            26 wfrct2
     7 mdl              17 score           27 mdl len
     8 mdl from         18 E-value         28 seq len
     9 mdl to           19 inc             29 description of target
    10 seq from         20 olp

For *cmscan* the CM is the target and the input sequence is the query, so
``target name``/``accession`` are the Rfam family ID and accession, and
``query name`` is the contig.  (This is the opposite of ``cmsearch``.)

Note that ``mdl len`` and ``seq len`` (27/28) are present in Infernal 1.1.5 but
not in the column layout given in some older documentation, so the description
starts at field 29.  Verified against real cmscan 1.1.5 output.
"""

import os

import polars as pl
from loguru import logger

# Column indices into a --fmt 2 cmscan tblout row.
_IDX_TARGET_NAME = 1
_IDX_TARGET_ACC = 2
_IDX_QUERY_NAME = 3
_IDX_CLAN = 5
_IDX_MDL_FROM = 7
_IDX_MDL_TO = 8
_IDX_SEQ_FROM = 9
_IDX_SEQ_TO = 10
_IDX_STRAND = 11
_IDX_TRUNC = 12
_IDX_GC = 14
_IDX_SCORE = 16
_IDX_EVALUE = 17
_IDX_INC = 18
_IDX_OLP = 19
# fields 20-25 are the overlap detail columns, 26 is mdl len and 27 is seq len;
# everything from field 28 on is the free-text description
_N_FIXED_FIELDS = 28

# tRNA and tmRNA are already annotated by tRNAscan-SE and ARAGORN respectively.
# Reporting them again from Rfam produces duplicate features in the GFF, and the
# specialised tools are more sensitive for phage sequence (Rfam's RF00023 misses
# phage tmRNAs that ARAGORN finds).
TRNA_TMRNA_ACCESSIONS = {"RF00005", "RF00023"}

NCRNA_TSV_COLUMNS = [
    "contig",
    "locus_tag",
    "start",
    "stop",
    "strand",
    "rfam_acc",
    "rfam_id",
    "type",
    "description",
    "clan",
    "bitscore",
    "evalue",
    "gc",
    "trunc",
    "mdl_from",
    "mdl_to",
]


# Rfam's free-text type strings (e.g. 'Cis-reg; riboswitch;') mapped onto the
# INSDC/GenBank ncRNA_class controlled vocabulary.  Checked in order, first
# match wins, so more specific terms must come first.  Anything unmatched
# becomes "other", which is a permitted INSDC value.
# NB: needles must be lowercase - they are matched against a lowercased type.
_NCRNA_CLASS_RULES = [
    ("snorna", "snoRNA"),
    ("snrna", "snRNA"),
    ("mirna", "miRNA"),
    ("ribozyme", "ribozyme"),
    ("antisense", "antisense_RNA"),
    ("tmrna", "tmRNA"),
    ("srp", "SRP_RNA"),
    ("rnasep", "RNase_P_RNA"),
    ("telomerase", "telomerase_RNA"),
    ("intron", "autocatalytically_spliced_intron"),
    ("lncrna", "lncRNA"),
    ("srna", "ncRNA"),
    ("vault", "vault_RNA"),
    ("y_rna", "Y_RNA"),
]


def rfam_type_to_ncrna_class(rfam_type):
    """Maps an Rfam type string onto an INSDC ncRNA_class value.

    NCBI requires ncRNA_class on every ncRNA feature, and only accepts terms
    from a fixed vocabulary.  Rfam's own type strings are free text and much
    broader (riboswitches, leaders, thermoregulators, frameshift elements),
    none of which have an ncRNA_class - those legitimately fall through to
    "other", with the specific Rfam family retained in the note attribute.
    """
    if not rfam_type:
        return "other"

    lowered = rfam_type.lower()
    for needle, ncrna_class in _NCRNA_CLASS_RULES:
        if needle in lowered:
            return ncrna_class

    return "other"


def load_rfam_metadata(db_dir):
    """Reads Rfam_metadata.tsv from the pharokka database directory.

    :param db_dir: pharokka database directory
    :return: dict keyed by Rfam accession (e.g. 'RF00023')
    """
    path = os.path.join(db_dir, "Rfam_metadata.tsv")

    metadata = {}
    df = pl.read_csv(path, separator="\t", has_header=True)
    for row in df.iter_rows(named=True):
        metadata[row["rfam_acc"]] = {
            "rfam_id": row.get("rfam_id", ""),
            "type": row.get("type", ""),
            "description": row.get("description", ""),
            "clan_acc": row.get("clan_acc", ""),
        }

    return metadata


def _empty_ncrna_df():
    """An empty dataframe with the ncRNA schema, so downstream code is uniform."""
    return pl.DataFrame(
        {
            "contig": pl.Series([], dtype=pl.Utf8),
            "locus_tag": pl.Series([], dtype=pl.Utf8),
            "start": pl.Series([], dtype=pl.Int64),
            "stop": pl.Series([], dtype=pl.Int64),
            "strand": pl.Series([], dtype=pl.Utf8),
            "rfam_acc": pl.Series([], dtype=pl.Utf8),
            "rfam_id": pl.Series([], dtype=pl.Utf8),
            "type": pl.Series([], dtype=pl.Utf8),
            "description": pl.Series([], dtype=pl.Utf8),
            "clan": pl.Series([], dtype=pl.Utf8),
            "bitscore": pl.Series([], dtype=pl.Float64),
            "evalue": pl.Series([], dtype=pl.Float64),
            "gc": pl.Series([], dtype=pl.Float64),
            "trunc": pl.Series([], dtype=pl.Utf8),
            "mdl_from": pl.Series([], dtype=pl.Int64),
            "mdl_to": pl.Series([], dtype=pl.Int64),
        }
    )


def parse_cmscan_tblout(tblout_path, metadata=None, keep_trna=False):
    """Parses a ``cmscan --fmt 2`` tblout into a tidy polars dataframe.

    Three filters are applied, in order:

    1. ``inc == '!'`` - the hit met the family's GA gathering threshold.
    2. ``olp != '='`` - drop hits marked as overlapping a higher-scoring hit
       from the same clan.  This is the clan competition step; without it,
       clans such as the riboswitches emit piles of redundant overlapping
       calls for the same locus.
    3. unless ``keep_trna``, drop tRNA/tmRNA families already covered by
       tRNAscan-SE and ARAGORN.

    :param tblout_path: path to the cmscan tblout
    :param metadata: dict from load_rfam_metadata(), or None
    :param keep_trna: keep RF00005/RF00023 hits
    :return: polars dataframe with NCRNA_TSV_COLUMNS (locus_tag left empty)
    """
    metadata = metadata or {}

    if not os.path.isfile(tblout_path):
        logger.warning(f"cmscan output {tblout_path} not found - no ncRNAs reported.")
        return _empty_ncrna_df()

    records = []
    with open(tblout_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split(None, _N_FIXED_FIELDS)
            if len(fields) <= _IDX_OLP:
                logger.warning(f"Skipping malformed cmscan line: {line[:80]}")
                continue

            if fields[_IDX_INC] != "!":
                continue
            if fields[_IDX_OLP] == "=":
                continue

            rfam_acc = fields[_IDX_TARGET_ACC]
            if not keep_trna and rfam_acc in TRNA_TMRNA_ACCESSIONS:
                continue

            # on the minus strand cmscan reports seq from > seq to
            seq_from = int(fields[_IDX_SEQ_FROM])
            seq_to = int(fields[_IDX_SEQ_TO])
            start, stop = min(seq_from, seq_to), max(seq_from, seq_to)

            clan = fields[_IDX_CLAN]
            info = metadata.get(rfam_acc, {})
            description = (
                fields[_N_FIXED_FIELDS] if len(fields) > _N_FIXED_FIELDS else ""
            )

            records.append(
                {
                    "contig": fields[_IDX_QUERY_NAME],
                    "locus_tag": "",
                    "start": start,
                    "stop": stop,
                    "strand": fields[_IDX_STRAND],
                    "rfam_acc": rfam_acc,
                    "rfam_id": info.get("rfam_id") or fields[_IDX_TARGET_NAME],
                    "type": info.get("type", ""),
                    "description": info.get("description") or description,
                    "clan": "" if clan == "-" else clan,
                    "bitscore": float(fields[_IDX_SCORE]),
                    "evalue": float(fields[_IDX_EVALUE]),
                    "gc": float(fields[_IDX_GC]),
                    "trunc": fields[_IDX_TRUNC],
                    "mdl_from": int(fields[_IDX_MDL_FROM]),
                    "mdl_to": int(fields[_IDX_MDL_TO]),
                }
            )

    if not records:
        return _empty_ncrna_df()

    df = pl.DataFrame(records, schema=_empty_ncrna_df().schema)
    return df.sort(["contig", "start"])


def add_locus_tags(df, locustag, contig_count):
    """Assigns pharokka locus tags to ncRNA rows.

    Follows the tRNA convention in post_processing.create_gff(): per-contig
    numbering when there are multiple contigs, otherwise a single run of
    numbers prefixed with the locustag.
    """
    if df.height == 0:
        return df

    if contig_count > 1:
        df = df.with_columns(
            (
                pl.col("contig")
                + pl.lit("_ncRNA_")
                + (pl.col("contig").cum_count().over("contig")).cast(pl.Utf8)
            ).alias("locus_tag")
        )
    else:
        df = df.with_columns(
            (
                pl.lit(f"{locustag}_ncRNA_")
                + (pl.int_range(1, pl.len() + 1)).cast(pl.Utf8)
            ).alias("locus_tag")
        )

    return df


def write_ncrna_tsv(df, out_dir, prefix):
    """Writes the {prefix}_ncrna.tsv output."""
    out_path = os.path.join(out_dir, f"{prefix}_ncrna.tsv")
    df.select(NCRNA_TSV_COLUMNS).write_csv(out_path, separator="\t")
    logger.info(f"{df.height} ncRNA(s) written to {out_path}")
