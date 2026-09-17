"""Unit tests for pure helpers in ``pharokka.post_processing``.

These tests target functions that don't need a database or a full pipeline
run.  Heavy ``Pharok`` class methods (``create_gff``, ``create_tbl``, etc.)
remain covered by the end-to-end integration tests in ``test_overall.py``.

Coverage focus:
- ``parse_attributes_column``      — vectorised attribute parser (fix #1)
- ``extract_anticodon_positions``  — dict-keyed lookup (bug fix #6)
- ``get_codon_from_anticodon``     — anticodon → codon RNA reverse-complement
- ``process_pyhmmer_results``      — pyhmmer hit → DataFrame columns
- ``process_custom_pyhmmer_results`` — same shape, custom HMM db
- ``process_vfdb_results``         — issue #410 bracket cleanup (fix #10)
- ``create_mmseqs_tophits``        — lowest-eVal-wins, empty input
- ``Pharok.parse_aragorn``        — out-of-contig ARAGORN spans (issue #444)
- ``is_trna_empty`` / ``is_file_empty`` / ``get_crispr_count`` / ``check_and_create_directory``
                                   — small filesystem helpers
"""

import collections
import warnings

import polars as pl
import pytest
from BCBio import GFF
from Bio import BiopythonParserWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pharokka.post_processing import (
    Pharok,
    check_and_create_directory,
    create_mmseqs_tophits,
    extract_anticodon_positions,
    get_codon_from_anticodon,
    get_crispr_count,
    is_file_empty,
    is_trna_empty,
    parse_attributes_column,
    process_custom_pyhmmer_results,
    process_pyhmmer_results,
    process_vfdb_results,
    read_feature_gff,
)

GFF_COLS = [
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

# ---------------------------------------------------------------------------
# parse_attributes_column
# ---------------------------------------------------------------------------


class TestParseAttributesColumn:
    """Vectorised GFF attribute parser — fix #1."""

    def test_basic_homogeneous(self):
        df = pl.DataFrame(
            {
                "contig": ["c1", "c1", "c2"],
                "attributes": [
                    "ID=a;phrog=1;product=foo",
                    "ID=b;phrog=2;product=bar",
                    "ID=c;phrog=3;product=baz",
                ],
            }
        )
        result = parse_attributes_column(df)
        # New columns appended in first-appearance order
        assert result.columns == ["contig", "attributes", "ID", "phrog", "product"]
        assert result["ID"].to_list() == ["a", "b", "c"]
        assert result["phrog"].to_list() == ["1", "2", "3"]
        assert result["product"].to_list() == ["foo", "bar", "baz"]

    def test_heterogeneous_keys_have_nulls(self):
        df = pl.DataFrame(
            {
                "attributes": [
                    "ID=a;phrog=1;partial=10",
                    "ID=b;phrog=2",
                    "ID=c;phrog=3;partial=01",
                ],
            }
        )
        result = parse_attributes_column(df)
        assert result.columns == ["attributes", "ID", "phrog", "partial"]
        assert result["partial"].to_list() == ["10", None, "01"]

    def test_column_order_preserves_first_appearance(self):
        """A key that first appears in row 1 (not 0) should still sit after
        keys from row 0."""
        df = pl.DataFrame(
            {
                "attributes": [
                    "ID=a;phrog=1",
                    "ID=b;phrog=2;partial=01;extra=foo",
                    "ID=c;phrog=3",
                ],
            }
        )
        result = parse_attributes_column(df)
        # extra appears AFTER partial because that's the order in row 1.
        assert result.columns == ["attributes", "ID", "phrog", "partial", "extra"]

    def test_empty_df_returns_unchanged(self):
        df = pl.DataFrame({"attributes": pl.Series([], dtype=pl.Utf8)})
        result = parse_attributes_column(df)
        assert result.height == 0
        assert result.columns == ["attributes"]

    def test_missing_attributes_column_returns_unchanged(self):
        df = pl.DataFrame({"contig": ["c1"]})
        result = parse_attributes_column(df)
        assert result.columns == ["contig"]


# ---------------------------------------------------------------------------
# extract_anticodon_positions  (regression test for bug fix #6)
# ---------------------------------------------------------------------------


class TestExtractAnticodonPositions:
    """Regression tests for the v1.9.1 → v1.10.0 anticodon-position bug fix.

    In v1.9.1 this returned a positional list, matched against ``trna_df``
    by row index in ``create_gff``.  The GFF is sorted by genomic position;
    the ``.sec`` file is sorted by tRNA ID — so positions were silently
    attached to the wrong tRNA whenever a contig had >1 tRNA.  The fix
    returns a dict keyed by tRNAscan ID.
    """

    SAMPLE_SEC = (
        "MW460250_1.trna1 (138642-138569)\tLength: 74 bp\n"
        "Type: Asp\tAnticodon: GTC at 35-37 (138608-138606)\tScore: 58.8\n"
        "Seq: GGCTCATTGGTGTAACTGGT\n"
        "Str: >>>>>>>..>>>>........\n"
        "\n"
        "MW460250_1.trna2 (138562-138490)\tLength: 73 bp\n"
        "Type: Phe\tAnticodon: GAA at 34-36 (138529-138527)\tScore: 67.8\n"
        "Seq: GGTTTCTTAGCTCAGATGGT\n"
        "Str: >>>>>>>..>>>>........\n"
        "\n"
        "MW460250_1.trna3 (115265-115194)\tLength: 72 bp\n"
        "Type: Met\tAnticodon: CAT at 33-35 (115233-115231)\tScore: 61.2\n"
        "Seq: GGACTCTTAGCTTAAAGGTA\n"
        "Str: >>>>>>>..>>>>.......\n"
    )

    def test_returns_dict_keyed_by_trnascan_id(self, tmp_path):
        (tmp_path / "trnascan_out.sec").write_text(self.SAMPLE_SEC)
        positions = extract_anticodon_positions(str(tmp_path))
        assert positions == {
            "MW460250_1.trna1": (138608, 138606),
            "MW460250_1.trna2": (138529, 138527),
            "MW460250_1.trna3": (115233, 115231),
        }

    def test_empty_sec_file_returns_empty_dict(self, tmp_path):
        (tmp_path / "trnascan_out.sec").write_text("")
        assert extract_anticodon_positions(str(tmp_path)) == {}

    def test_single_trna(self, tmp_path):
        sec = (
            "NC_051700.1.trna1 (32106-32041)\tLength: 66 bp\n"
            "Type: Lys\tAnticodon: CTT at 28-30 (32079-32077)\tScore: 17.1\n"
            "Seq: TGGTatGCAGCCTATCaACCG\n"
        )
        (tmp_path / "trnascan_out.sec").write_text(sec)
        positions = extract_anticodon_positions(str(tmp_path))
        assert positions == {"NC_051700.1.trna1": (32079, 32077)}

    def test_ignores_non_header_non_type_lines(self, tmp_path):
        """Blank lines, Seq:, and Str: lines must not be parsed as headers."""
        (tmp_path / "trnascan_out.sec").write_text(self.SAMPLE_SEC)
        positions = extract_anticodon_positions(str(tmp_path))
        # Same expected as above — confirms the Seq:/Str: lines weren't picked up.
        assert len(positions) == 3


# ---------------------------------------------------------------------------
# get_codon_from_anticodon
# ---------------------------------------------------------------------------


class TestGetCodonFromAnticodon:
    @pytest.mark.parametrize(
        ("anticodon", "expected_codon"),
        [
            ("CAT", "AUG"),  # Met
            ("CTT", "AAG"),  # Lys
            ("GAA", "UUC"),  # Phe
            ("GTC", "GAC"),  # Asp
            ("cat", "AUG"),  # lowercase input is uppercased
        ],
    )
    def test_known_codons(self, anticodon, expected_codon):
        assert get_codon_from_anticodon(anticodon) == expected_codon

    def test_non_string_raises_type_error(self):
        with pytest.raises(TypeError):
            get_codon_from_anticodon(None)

    def test_empty_string_raises_value_error(self):
        with pytest.raises(ValueError):
            get_codon_from_anticodon("")

    def test_whitespace_only_raises_value_error(self):
        with pytest.raises(ValueError):
            get_codon_from_anticodon("   ")


# ---------------------------------------------------------------------------
# process_pyhmmer_results
# ---------------------------------------------------------------------------

# Match the namedtuple definition at top of post_processing.py
Result = collections.namedtuple("Result", ["protein", "phrog", "bitscore", "evalue"])
CustomResult = collections.namedtuple(
    "Result", ["protein", "custom_hmm_id", "bitscore", "evalue"]
)


class TestProcessPyhmmerResults:
    def _base_df(self):
        return pl.DataFrame(
            {
                "gene": ["gene_001 100_200", "gene_002 300_400", "gene_003 500_600"],
                "other_col": ["x", "y", "z"],
            }
        )

    def test_attaches_hits_for_matching_proteins(self):
        df = self._base_df()
        # Gene IDs in pyhmmer_results_dict use the pre-space portion of "gene"
        hits = {
            "gene_001": Result("gene_001", "phrog_42", 100.5, 1e-30),
            "gene_003": Result("gene_003", "phrog_7", 50.25, 1e-10),
        }
        result = process_pyhmmer_results(df, hits)
        assert "pyhmmer_phrog" in result.columns
        assert "pyhmmer_bitscore" in result.columns
        assert "pyhmmer_evalue" in result.columns
        assert result["pyhmmer_phrog"].to_list() == [
            "phrog_42",
            "No_PHROGs_HMM",
            "phrog_7",
        ]
        # bitscore rounded to 6 decimal places, str()-ified
        assert result["pyhmmer_bitscore"].to_list() == [
            "100.5",
            "No_PHROGs_HMM",
            "50.25",
        ]
        assert result["pyhmmer_evalue"].to_list() == ["1e-30", "No_PHROGs_HMM", "1e-10"]

    def test_empty_hits_dict_yields_all_no_phrogs(self):
        df = self._base_df()
        result = process_pyhmmer_results(df, {})
        assert result["pyhmmer_phrog"].to_list() == ["No_PHROGs_HMM"] * 3
        assert result["pyhmmer_bitscore"].to_list() == ["No_PHROGs_HMM"] * 3
        assert result["pyhmmer_evalue"].to_list() == ["No_PHROGs_HMM"] * 3

    def test_gene_name_split_strips_post_space_portion(self):
        """The protein name used for lookup is everything before the first space."""
        df = pl.DataFrame({"gene": ["abc def ghi"]})
        result = process_pyhmmer_results(
            df, {"abc": Result("abc", "phrog_1", 1.0, 1e-5)}
        )
        assert result["pyhmmer_phrog"].to_list() == ["phrog_1"]


class TestProcessCustomPyhmmerResults:
    def test_attaches_custom_hits(self):
        df = pl.DataFrame({"gene": ["g1 100_200", "g2 300_400"]})
        hits = {"g1": CustomResult("g1", "custom_xyz", 99.9, 1e-25)}
        result = process_custom_pyhmmer_results(df, hits)
        assert result["custom_hmm_id"].to_list() == ["custom_xyz", "No_custom_HMM"]
        assert result["custom_hmm_bitscore"].to_list() == ["99.9", "No_custom_HMM"]
        assert result["custom_hmm_evalue"].to_list() == ["1e-25", "No_custom_HMM"]


# ---------------------------------------------------------------------------
# process_vfdb_results  (fix #10: issue #410 bracket cleanup)
# ---------------------------------------------------------------------------


class TestProcessVfdbResults:
    """Confirm the 8 chained bracket-cleanup substitutions produce the same
    output as v1.9.1's 8 sequential ones."""

    @pytest.fixture
    def merged_df(self):
        # process_vfdb_results joins by "gene" onto merged_df
        return pl.DataFrame({"gene": ["gene_001", "gene_002", "gene_003"]})

    def test_empty_vfdb_file_no_results(self, tmp_path, merged_df):
        """When the VFDB results file is missing/empty, tophits is empty
        but the function should still return without raising."""
        # touch_file is called inside the function for the missing case
        merged_df_out, tophits = process_vfdb_results(
            str(tmp_path), merged_df, proteins_flag=False, reverse_mmseqs2=False
        )
        assert tophits.height == 0

    def test_bracket_substitutions_applied(self, tmp_path, merged_df):
        """Each of the 8 patterns from issue #410 / GenBank parsing fixes
        gets cleaned up."""
        # Build a synthetic vfdb_results.tsv with one row per pattern.
        # Column order (forward mode): vfdb_hit, gene, alnScore, seqIdentity, eVal,
        #                              qStart, qEnd, qLen, tStart, tEnd, tLen
        rows = [
            (
                "VFG001 L-allo-isoleucine:holo-[CmaA peptidyl-carrier protein] foo",
                "gene_001",
                "100",
                "0.9",
                "1e-30",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
            (
                "VFG002 UDP-3-O-[3-hydroxymyristoyl] synthase",
                "gene_002",
                "200",
                "0.8",
                "1e-25",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
            (
                "VFG003 biotin--[acetyl-CoA-carboxylase] protein",
                "gene_003",
                "300",
                "0.7",
                "1e-20",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
        ]
        vfdb_file = tmp_path / "vfdb_results.tsv"
        with open(vfdb_file, "w") as f:
            for row in rows:
                f.write("\t".join(row) + "\n")

        _, tophits = process_vfdb_results(
            str(tmp_path), merged_df, proteins_flag=False, reverse_mmseqs2=False
        )
        hits = tophits["vfdb_hit"].to_list()
        # bracket characters should be stripped from these specific patterns
        assert any(
            "L-allo-isoleucine:holo-CmaA peptidyl-carrier protein" in h for h in hits
        )
        assert any("UDP-3-O-3-hydroxymyristoyl synthase" in h for h in hits)
        assert any("biotin--acetyl-CoA-carboxylase protein" in h for h in hits)
        # No leftover square brackets on the cleaned-up patterns
        assert not any("[CmaA" in h for h in hits)
        assert not any("UDP-3-O-[3-hydroxymyristoyl] synthase" in h for h in hits)

    def test_lowest_evalue_wins_per_gene(self, tmp_path, merged_df):
        """Two hits for the same gene → the one with the smallest eVal wins."""
        rows = [
            (
                "VFG_LOSER",
                "gene_001",
                "100",
                "0.9",
                "1e-10",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
            (
                "VFG_WINNER",
                "gene_001",
                "200",
                "0.8",
                "1e-30",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
        ]
        vfdb_file = tmp_path / "vfdb_results.tsv"
        with open(vfdb_file, "w") as f:
            for row in rows:
                f.write("\t".join(row) + "\n")

        _, tophits = process_vfdb_results(
            str(tmp_path), merged_df, proteins_flag=False, reverse_mmseqs2=False
        )
        gene_001_rows = tophits.filter(pl.col("gene") == "gene_001")
        assert gene_001_rows.height == 1
        assert "VFG_WINNER" in gene_001_rows["vfdb_hit"].to_list()[0]


# ---------------------------------------------------------------------------
# create_mmseqs_tophits
# ---------------------------------------------------------------------------


class TestCreateMmseqsTophits:
    def test_empty_mmseqs_file_returns_empty_df(self, tmp_path):
        """When the mmseqs file is empty, polars raises ``NoDataError`` —
        the function catches it and returns an empty DataFrame."""
        (tmp_path / "mmseqs_results.tsv").write_text("")
        tophits = create_mmseqs_tophits(str(tmp_path), reverse_mmseqs=False)
        assert tophits.height == 0

    def test_picks_lowest_evalue_per_gene(self, tmp_path):
        """Two hits for the same gene → the lower-eVal hit wins."""
        # Forward mode columns: phrog, gene, alnScore, seqIdentity, eVal,
        #                      qStart, qEnd, qLen, tStart, tEnd, tLen
        rows = [
            (
                "phrog_loser",
                "gene_001",
                "100",
                "0.9",
                "1e-10",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
            (
                "phrog_winner",
                "gene_001",
                "200",
                "0.8",
                "1e-30",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
            (
                "phrog_other",
                "gene_002",
                "50",
                "0.7",
                "1e-5",
                "1",
                "100",
                "100",
                "1",
                "100",
                "100",
            ),
        ]
        with open(tmp_path / "mmseqs_results.tsv", "w") as f:
            for row in rows:
                f.write("\t".join(row) + "\n")

        tophits = create_mmseqs_tophits(str(tmp_path), reverse_mmseqs=False)
        assert tophits.height == 2  # one row per unique gene
        g1 = tophits.filter(pl.col("gene") == "gene_001")
        assert g1["mmseqs_phrog"].to_list() == ["phrog_winner"]


# ---------------------------------------------------------------------------
# Small filesystem helpers
# ---------------------------------------------------------------------------


class TestFilesystemHelpers:
    def test_is_file_empty_true_for_empty_file(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("")
        assert is_file_empty(str(empty)) is True

    def test_is_file_empty_false_for_nonempty_file(self, tmp_path):
        nonempty = tmp_path / "data.txt"
        nonempty.write_text("hello")
        assert is_file_empty(str(nonempty)) is False

    def test_is_trna_empty_true_when_gff_empty(self, tmp_path):
        (tmp_path / "trnascan_out.gff").write_text("")
        assert is_trna_empty(str(tmp_path)) is True

    def test_is_trna_empty_false_when_gff_has_content(self, tmp_path):
        (tmp_path / "trnascan_out.gff").write_text("##gff-version 3\n")
        assert is_trna_empty(str(tmp_path)) is False

    def test_get_crispr_count_zero(self, tmp_path):
        (tmp_path / "test_minced.gff").write_text("##gff-version 3\n# comment line\n")
        assert get_crispr_count(str(tmp_path), "test") == 0

    def test_get_crispr_count_counts_non_comment_lines(self, tmp_path):
        content = (
            "##gff-version 3\n"
            "# more comment\n"
            "contig1\tminced\trepeat_region\t100\t200\t.\t+\t.\tID=CRISPR1\n"
            "contig1\tminced\trepeat_region\t300\t400\t.\t+\t.\tID=CRISPR2\n"
        )
        (tmp_path / "test_minced.gff").write_text(content)
        assert get_crispr_count(str(tmp_path), "test") == 2

    def test_check_and_create_directory_creates_if_missing(self, tmp_path):
        new_dir = tmp_path / "newdir"
        assert not new_dir.exists()
        check_and_create_directory(str(new_dir))
        assert new_dir.is_dir()

    def test_check_and_create_directory_noop_if_exists(self, tmp_path):
        existing = tmp_path / "existing"
        existing.mkdir()
        # Should not raise
        check_and_create_directory(str(existing))
        assert existing.is_dir()


class TestReadFeatureGff:
    """Regression tests for read_feature_gff.

    A phage with no tRNAs (or CRISPRs / tmRNAs) leaves a 0-byte feature GFF.
    Passing that to polars raised NoDataError on some platforms and
    OSError ("Exec format error") on others, crashing post-processing.
    read_feature_gff must return an empty frame instead of raising.
    """

    def test_zero_byte_file_returns_empty_frame(self, tmp_path):
        gff = tmp_path / "trnascan_out.gff"
        gff.write_text("")
        df = read_feature_gff(str(gff), GFF_COLS)
        assert df.height == 0
        assert df.columns == GFF_COLS

    def test_header_only_file_returns_empty_frame(self, tmp_path):
        gff = tmp_path / "trnascan_out.gff"
        gff.write_text("##gff-version 3\n")
        df = read_feature_gff(str(gff), GFF_COLS)
        assert df.height == 0
        assert df.columns == GFF_COLS

    def test_missing_file_returns_empty_frame(self, tmp_path):
        df = read_feature_gff(str(tmp_path / "does_not_exist.gff"), GFF_COLS)
        assert df.height == 0
        assert df.columns == GFF_COLS

    def test_reads_feature_rows_and_skips_comments(self, tmp_path):
        gff = tmp_path / "trnascan_out.gff"
        gff.write_text(
            "##gff-version 3\n"
            "ctg1\ttRNAscan-SE\ttRNA\t10\t80\t.\t+\t.\tID=ctg1.trna1\n"
            "ctg1\ttRNAscan-SE\ttRNA\t120\t190\t.\t-\t.\tID=ctg1.trna2\n"
        )
        df = read_feature_gff(str(gff), GFF_COLS)
        assert df.height == 2
        assert df["Region"].to_list() == ["tRNA", "tRNA"]
        assert df["contig"][0] == "ctg1"


# ---------------------------------------------------------------------------
# parse_aragorn — out-of-contig coordinates (issue #444, phold #141)
# ---------------------------------------------------------------------------


def _write_aragorn(tmp_path, prefix, body):
    (tmp_path / f"{prefix}_aragorn.txt").write_text(body)


def _aragorn_gff_rows(tmp_path, prefix):
    text = (tmp_path / f"{prefix}_aragorn.gff").read_text().strip()
    return [line.split("\t") for line in text.splitlines() if line.strip()]


def _run_parse_aragorn(tmp_path, prefix, body, contigs, lengths, meta_mode=False):
    _write_aragorn(tmp_path, prefix, body)
    pharok = Pharok(
        out_dir=str(tmp_path),
        prefix=prefix,
        aragorn_version="1.2.41",
        meta_mode=meta_mode,
        length_df=pl.DataFrame({"contig": contigs, "length": lengths}),
    )
    pharok.parse_aragorn()
    return pharok


class TestParseAragornCoordinates:
    """ARAGORN spans that run off the end of a contig — issue #444.

    Run linearly (``-l``), ARAGORN reports a gene whose model overhangs the 5'
    end of a sequence with a non-positive start, stepping over zero so the base
    before position 1 is ``-1`` (``position()`` in aragorn.c).  Written through
    to the GenBank unchanged, a location like ``-63..456`` makes BioPython
    *warn* and set ``SeqFeature.location`` to ``None``, and every downstream
    tool that sorts on ``feature.location.start`` then dies on the ``None``
    (phold #141).  parse_aragorn must clamp such spans onto the contig.
    """

    def test_negative_start_is_clamped_to_one(self, tmp_path):
        _run_parse_aragorn(
            tmp_path,
            "neg",
            ">c1\n1 gene found\n1   tmRNA  [-63,456]\t157,258\tASARS*\n",
            ["c1"],
            [5000],
        )
        (row,) = _aragorn_gff_rows(tmp_path, "neg")
        assert row[3:5] == ["1", "456"]

    def test_stop_past_contig_end_is_clamped(self, tmp_path):
        _run_parse_aragorn(
            tmp_path,
            "over",
            ">c1\n1 gene found\n1   tmRNA  [4800,5200]\t157,258\tASARS*\n",
            ["c1"],
            [5000],
        )
        (row,) = _aragorn_gff_rows(tmp_path, "over")
        assert row[3:5] == ["4800", "5000"]

    def test_in_range_span_is_untouched(self, tmp_path):
        pharok = _run_parse_aragorn(
            tmp_path,
            "ok",
            ">c1\n1 gene found\n1   tmRNA  [100,454]\t157,258\tASARS*\n",
            ["c1"],
            [5000],
        )
        (row,) = _aragorn_gff_rows(tmp_path, "ok")
        assert row[3:5] == ["100", "454"]
        assert pharok.tmrna_flag is True

    def test_complementary_strand_negative_start_is_clamped(self, tmp_path):
        """The 'c' strand marker is stripped before the coordinate is read."""
        _run_parse_aragorn(
            tmp_path,
            "comp",
            ">c1\n1 gene found\n1   tmRNA  c[-10,300]\t157,258\tASARS*\n",
            ["c1"],
            [5000],
        )
        (row,) = _aragorn_gff_rows(tmp_path, "comp")
        assert row[3:5] == ["1", "300"]

    def test_span_entirely_off_contig_is_dropped(self, tmp_path):
        """Nothing survives, so tmrna_flag must go back to False.

        Left True, create_gff() would go on to read the now-empty
        {prefix}_aragorn.gff and raise.
        """
        pharok = _run_parse_aragorn(
            tmp_path,
            "off",
            ">c1\n1 gene found\n1   tmRNA  [-500,-100]\t157,258\tASARS*\n",
            ["c1"],
            [5000],
        )
        assert _aragorn_gff_rows(tmp_path, "off") == []
        assert pharok.tmrna_flag is False

    def test_multi_contig_clamps_against_the_right_contig(self, tmp_path):
        """Each span is clamped against its own contig's length, not the first."""
        pharok = _run_parse_aragorn(
            tmp_path,
            "multi",
            ">c1\n1 gene found\n1   tmRNA  [-63,456]\t157,258\tASARS*\n"
            ">c2\n1 gene found\n1   tmRNA  [100,454]\t157,258\tASARS*\n"
            ">c3\n1 gene found\n1   tmRNA  [400,900]\t157,258\tASARS*\n"
            ">c4\n0 genes found\n"
            ">end \t4 sequences 3 tmRNA genes, nothing found in 1 sequences,"
            " (75.00% sensitivity)\n",
            ["c1", "c2", "c3", "c4"],
            [5000, 5000, 600, 5000],
            meta_mode=True,
        )
        rows = _aragorn_gff_rows(tmp_path, "multi")
        assert [(r[0], r[3], r[4]) for r in rows] == [
            ("c1", "1", "456"),
            ("c2", "100", "454"),
            ("c3", "400", "600"),
        ]
        assert pharok.tmrna_flag is True

    def test_clamped_span_survives_the_round_trip_to_genbank(self, tmp_path):
        """The end-to-end chain from phold #141, pinned at the pharokka end.

        Before the fix the GFF carried ``-63  456``, BioPython wrote the
        location ``-63..456`` and re-read it as ``None``.
        """
        _run_parse_aragorn(
            tmp_path,
            "gbk",
            ">c1\n1 gene found\n1   tmRNA  [-63,456]\t157,258\tASARS*\n",
            ["c1"],
            [5000],
        )
        (row,) = _aragorn_gff_rows(tmp_path, "gbk")

        gff = tmp_path / "mini.gff"
        gff.write_text(
            "##gff-version 3\n"
            f"c1\t{row[1]}\ttmRNA\t{row[3]}\t{row[4]}\t.\t+\t.\tID=c1_tmRNA_0001\n"
        )
        contig = SeqRecord(Seq("A" * 5000), id="c1")
        (record,) = list(GFF.parse(str(gff), {"c1": contig}))
        record.annotations.update(
            molecule_type="DNA", topology="linear", data_file_division="PHG"
        )
        gbk = tmp_path / "mini.gbk"
        SeqIO.write(record, str(gbk), "genbank")

        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            (reparsed,) = list(SeqIO.parse(str(gbk), "genbank"))

        (feature,) = reparsed.features
        assert feature.location is not None
        # 0-based half-open, so GFF 1..456 -> [0:456]
        assert (int(feature.location.start), int(feature.location.end)) == (0, 456)


class TestParseAragornEndSummary:
    """ARAGORN's trailing ``>end`` summary line must not be read as a contig.

    On a multi-sequence run ARAGORN closes its batch output with
    ``>end \\t<N> sequences ...``.  That line only mentions "sensitivity" when
    at least one sequence came up empty, so a run where *every* contig carried
    a tmRNA used to fall through the header check, be treated as one more
    contig, and index ``lines[i + 1]`` off the end of the file.
    """

    def test_every_contig_with_a_tmrna_does_not_crash(self, tmp_path):
        pharok = _run_parse_aragorn(
            tmp_path,
            "allhit",
            ">c1\n1 gene found\n1   tmRNA  [100,454]\t157,258\tASARS*\n"
            ">c2\n1 gene found\n1   tmRNA  [200,554]\t157,258\tASARS*\n"
            ">end \t2 sequences 2 tmRNA genes\n",
            ["c1", "c2"],
            [5000, 5000],
            meta_mode=True,
        )
        rows = _aragorn_gff_rows(tmp_path, "allhit")
        assert [(r[0], r[3], r[4]) for r in rows] == [
            ("c1", "100", "454"),
            ("c2", "200", "554"),
        ]
        assert pharok.tmrna_flag is True

    def test_summary_with_sensitivity_still_skipped(self, tmp_path):
        """The wording ARAGORN uses when some contig came up empty."""
        pharok = _run_parse_aragorn(
            tmp_path,
            "somehit",
            ">c1\n1 gene found\n1   tmRNA  [100,454]\t157,258\tASARS*\n"
            ">c2\n0 genes found\n"
            ">end \t2 sequences 1 tmRNA genes, nothing found in 1 sequences,"
            " (50.00% sensitivity)\n",
            ["c1", "c2"],
            [5000, 5000],
            meta_mode=True,
        )
        rows = _aragorn_gff_rows(tmp_path, "somehit")
        assert [(r[0], r[3], r[4]) for r in rows] == [("c1", "100", "454")]
        assert pharok.tmrna_flag is True

    def test_contig_named_end_is_still_parsed(self, tmp_path):
        """Bounding by contig count, not by the ">end" text, so this is safe."""
        _run_parse_aragorn(
            tmp_path,
            "endname",
            ">c1\n0 genes found\n"
            ">end_of_assembly\n1 gene found\n1   tmRNA  [10,364]\t157,258\tASARS*\n"
            ">end \t2 sequences 1 tmRNA genes, nothing found in 1 sequences,"
            " (50.00% sensitivity)\n",
            ["c1", "end_of_assembly"],
            [5000, 5000],
            meta_mode=True,
        )
        rows = _aragorn_gff_rows(tmp_path, "endname")
        assert [(r[0], r[3], r[4]) for r in rows] == [("end_of_assembly", "10", "364")]
