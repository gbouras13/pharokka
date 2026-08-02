"""Tests for Rfam / Infernal ncRNA annotation."""

import os
import shutil
import subprocess

import polars as pl
import pytest

from pharokka.databases import RFAM_DB_NAMES, check_rfam_installation
from pharokka.post_processing import Pharok
from pharokka.rfam import (
    TRNA_TMRNA_ACCESSIONS,
    add_locus_tags,
    load_rfam_metadata,
    parse_cmscan_tblout,
    rfam_type_to_ncrna_class,
    write_ncrna_tsv,
)

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")
RFAM_DIR = os.path.join(TEST_DATA, "rfam")
NC_004617 = os.path.join(TEST_DATA, "overall", "VFDB_example", "NC_004617.fasta")
NC_051700 = os.path.join(TEST_DATA, "overall", "tmRNA_example", "NC_051700.fasta")

# A real cmscan --fmt 2 header, kept verbatim so the fixtures below stay
# aligned with the column layout the parser assumes.
TBLOUT_HEADER = """\
#idx target name          accession query name           accession clan name mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc olp anyidx afrct1 afrct2 winidx wfrct1 wfrct2 mdl len seq len description of target
#--- -------------------- --------- -------------------- --------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- --- ------ ------ ------ ------ ------ ------ ------- ------- ---------------------
"""


def _row(
    idx=1,
    name="SprX",
    acc="RF02672",
    query="NC_004617.1",
    clan="-",
    seq_from=39255,
    seq_to=39405,
    strand="+",
    score=195.7,
    evalue="1.2e-44",
    inc="!",
    olp="*",
    desc="Small pathogenicity island RNA X",
):
    return (
        f"{idx}    {name}    {acc}    {query}    -    {clan}    cm    1    151    "
        f"{seq_from}    {seq_to}    {strand}    no    1  0.37   0.1  {score}   "
        f"{evalue}   {inc}   {olp}   -      -      -      -      -      -     151   "
        f"42722 {desc}\n"
    )


def _write_tblout(tmp_path, rows):
    path = tmp_path / "test_cmscan.tblout"
    path.write_text(TBLOUT_HEADER + "".join(rows) + "#\n")
    return str(path)


class TestParseCmscanTblout:
    def test_parses_a_basic_hit(self, tmp_path):
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row()]))

        assert df.height == 1
        row = df.row(0, named=True)
        assert row["contig"] == "NC_004617.1"
        assert row["rfam_acc"] == "RF02672"
        assert row["rfam_id"] == "SprX"
        assert row["start"] == 39255
        assert row["stop"] == 39405
        assert row["strand"] == "+"
        assert row["bitscore"] == pytest.approx(195.7)
        assert row["evalue"] == pytest.approx(1.2e-44)

    def test_description_is_read_from_the_final_field(self, tmp_path):
        """Guards the field count: mdl len and seq len precede the description."""
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row()]))
        assert df.row(0, named=True)["description"] == (
            "Small pathogenicity island RNA X"
        )

    def test_below_gathering_threshold_hits_are_dropped(self, tmp_path):
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row(inc="?")]))
        assert df.height == 0

    def test_clan_competition_losers_are_dropped(self, tmp_path):
        """olp '=' marks a hit overlapping a better one in the same clan."""
        rows = [_row(idx=1, olp="^"), _row(idx=2, acc="RF00050", olp="=")]
        df = parse_cmscan_tblout(_write_tblout(tmp_path, rows))

        assert df.height == 1
        assert df.row(0, named=True)["rfam_acc"] == "RF02672"

    def test_minus_strand_coordinates_are_normalised(self, tmp_path):
        """cmscan reports seq from > seq to on the minus strand."""
        rows = [_row(seq_from=138562, seq_to=138490, strand="-")]
        df = parse_cmscan_tblout(_write_tblout(tmp_path, rows))

        row = df.row(0, named=True)
        assert row["start"] == 138490
        assert row["stop"] == 138562
        assert row["strand"] == "-"

    @pytest.mark.parametrize("accession", sorted(TRNA_TMRNA_ACCESSIONS))
    def test_trna_and_tmrna_dropped_by_default(self, tmp_path, accession):
        """tRNAscan-SE and ARAGORN already annotate these."""
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row(acc=accession)]))
        assert df.height == 0

    @pytest.mark.parametrize("accession", sorted(TRNA_TMRNA_ACCESSIONS))
    def test_trna_and_tmrna_kept_with_flag(self, tmp_path, accession):
        df = parse_cmscan_tblout(
            _write_tblout(tmp_path, [_row(acc=accession)]), keep_trna=True
        )
        assert df.height == 1

    def test_hits_are_sorted_by_position(self, tmp_path):
        rows = [
            _row(idx=1, seq_from=40870, seq_to=41011, acc="RF01828"),
            _row(idx=2, seq_from=36836, seq_to=37015, acc="RF01492"),
        ]
        df = parse_cmscan_tblout(_write_tblout(tmp_path, rows))
        assert df["start"].to_list() == [36836, 40870]

    def test_missing_file_returns_empty_frame(self, tmp_path):
        df = parse_cmscan_tblout(str(tmp_path / "nope.tblout"))
        assert df.height == 0
        assert "rfam_acc" in df.columns

    def test_empty_result_keeps_schema(self, tmp_path):
        df = parse_cmscan_tblout(_write_tblout(tmp_path, []))
        assert df.height == 0
        assert (
            df.columns == parse_cmscan_tblout(_write_tblout(tmp_path, [_row()])).columns
        )

    def test_metadata_overrides_target_name(self, tmp_path):
        metadata = {
            "RF02672": {
                "rfam_id": "SprX",
                "type": "Gene; sRNA;",
                "description": "Small pathogenicity island RNA X",
                "clan_acc": "",
            }
        }
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row()]), metadata=metadata)
        assert df.row(0, named=True)["type"] == "Gene; sRNA;"


class TestNcrnaClass:
    @pytest.mark.parametrize(
        "rfam_type,expected",
        [
            ("Gene; sRNA;", "ncRNA"),
            ("Gene; ribozyme;", "ribozyme"),
            ("Gene; antisense;", "antisense_RNA"),
            ("Gene; snRNA; snoRNA; CD-box;", "snoRNA"),
            ("Gene; snRNA; splicing;", "snRNA"),
            ("Gene; miRNA;", "miRNA"),
            ("Intron;", "autocatalytically_spliced_intron"),
            # riboswitches have no INSDC ncRNA_class
            ("Cis-reg; riboswitch;", "other"),
            ("Cis-reg; thermoregulator;", "other"),
            ("", "other"),
            (None, "other"),
        ],
    )
    def test_mapping(self, rfam_type, expected):
        assert rfam_type_to_ncrna_class(rfam_type) == expected


class TestLocusTags:
    def test_single_contig_uses_locustag_prefix(self, tmp_path):
        rows = [
            _row(idx=1, seq_from=100, seq_to=200),
            _row(idx=2, seq_from=300, seq_to=400),
        ]
        df = add_locus_tags(
            parse_cmscan_tblout(_write_tblout(tmp_path, rows)), "ABCDE", 1
        )
        assert df["locus_tag"].to_list() == ["ABCDE_ncRNA_1", "ABCDE_ncRNA_2"]

    def test_multi_contig_numbers_per_contig(self, tmp_path):
        rows = [
            _row(idx=1, query="contig1", seq_from=100, seq_to=200),
            _row(idx=2, query="contig1", seq_from=300, seq_to=400),
            _row(idx=3, query="contig2", seq_from=100, seq_to=200),
        ]
        df = add_locus_tags(
            parse_cmscan_tblout(_write_tblout(tmp_path, rows)), "ABCDE", 2
        )
        assert df["locus_tag"].to_list() == [
            "contig1_ncRNA_1",
            "contig1_ncRNA_2",
            "contig2_ncRNA_1",
        ]

    def test_empty_frame_is_left_alone(self, tmp_path):
        df = add_locus_tags(
            parse_cmscan_tblout(_write_tblout(tmp_path, [])), "ABCDE", 1
        )
        assert df.height == 0


class TestMetadataAndOutput:
    def test_load_metadata(self):
        metadata = load_rfam_metadata(RFAM_DIR)
        assert metadata["RF00023"]["rfam_id"] == "tmRNA"
        assert metadata["RF01828"]["type"] == "Gene; sRNA;"

    def test_write_ncrna_tsv(self, tmp_path):
        df = add_locus_tags(
            parse_cmscan_tblout(_write_tblout(tmp_path, [_row()])), "ABCDE", 1
        )
        write_ncrna_tsv(df, str(tmp_path), "test")

        out = tmp_path / "test_ncrna.tsv"
        lines = out.read_text().strip().split("\n")
        assert lines[0].split("\t")[:3] == ["contig", "locus_tag", "start"]
        assert len(lines) == 2

    def test_write_ncrna_tsv_when_empty(self, tmp_path):
        """A header-only file is still written, so the output set is stable."""
        write_ncrna_tsv(
            parse_cmscan_tblout(str(tmp_path / "no.tblout")), str(tmp_path), "t"
        )
        assert (
            (tmp_path / "t_ncrna.tsv")
            .read_text()
            .strip()
            .split("\n")[0]
            .startswith("contig")
        )


class TestRfamDatabaseCheck:
    """check_rfam_installation gates the default-on ncRNA path."""

    def test_passes_on_a_complete_database(self, tmp_path):
        for name in RFAM_DB_NAMES:
            (tmp_path / name).touch()
        assert check_rfam_installation(str(tmp_path)) is True

    def test_exits_when_a_file_is_missing(self, tmp_path):
        """pharokka's logger.error sink exits, so this never returns False."""
        for name in RFAM_DB_NAMES[:-1]:
            (tmp_path / name).touch()
        with pytest.raises(SystemExit):
            check_rfam_installation(str(tmp_path))

    def test_exits_on_an_empty_directory(self, tmp_path):
        with pytest.raises(SystemExit):
            check_rfam_installation(str(tmp_path))

    def test_flatfile_is_not_required(self):
        """cmscan reads the pressed .i1* files; shipping Rfam.cm would waste 329 MB."""
        assert "Rfam.cm" not in RFAM_DB_NAMES
        assert "Rfam.cm.i1m" in RFAM_DB_NAMES


class TestGffConstruction:
    """Pharok._build_ncrna_gff_df turns parsed hits into GFF3 rows."""

    def _pharok(self, ncrna_df):
        pharok = Pharok()
        pharok.ncrna_df = ncrna_df
        pharok.locustag = "ABCDE"
        pharok.infernal_version = "1.1.5"
        pharok.length_df = pl.DataFrame({"contig": ["NC_004617.1"], "length": [42722]})
        return pharok

    def test_gff_columns_and_values(self, tmp_path):
        metadata = {
            "RF02672": {
                "rfam_id": "SprX",
                "type": "Gene; sRNA;",
                "description": "Small pathogenicity island RNA X",
                "clan_acc": "",
            }
        }
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row()]), metadata=metadata)
        gff = self._pharok(df)._build_ncrna_gff_df()

        assert gff.columns == [
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

        row = gff.row(0, named=True)
        assert row["Region"] == "ncRNA"
        assert row["Method"] == "profile:Infernal:1.1.5"
        assert row["frame"] == "."
        assert row["attributes"] == (
            "ID=ABCDE_ncRNA_1;locus_tag=ABCDE_ncRNA_1;"
            "product=Small pathogenicity island RNA X;"
            "Dbxref=RFAM:RF02672;ncRNA_class=ncRNA;note=SprX"
        )

    def test_locus_tags_are_written_back_for_the_tsv(self, tmp_path):
        """The GFF and _ncrna.tsv must agree on locus tags."""
        df = parse_cmscan_tblout(_write_tblout(tmp_path, [_row()]))
        pharok = self._pharok(df)
        pharok._build_ncrna_gff_df()

        assert pharok.ncrna_df["locus_tag"].to_list() == ["ABCDE_ncRNA_1"]

    def test_riboswitch_gets_other_ncrna_class(self, tmp_path):
        metadata = {
            "RF00050": {
                "rfam_id": "FMN",
                "type": "Cis-reg; riboswitch;",
                "description": "FMN riboswitch",
                "clan_acc": "",
            }
        }
        df = parse_cmscan_tblout(
            _write_tblout(tmp_path, [_row(acc="RF00050")]), metadata=metadata
        )
        gff = self._pharok(df)._build_ncrna_gff_df()
        assert "ncRNA_class=other" in gff.row(0, named=True)["attributes"]


@pytest.fixture(scope="module")
def pressed_db(tmp_path_factory):
    """cmpress the committed 5-model subset into a temp database directory.

    The subset is committed unpressed (the .i1* files are large and binary);
    pressing takes well under a second for five models.
    """
    db = tmp_path_factory.mktemp("rfam_db")
    for name in ("Rfam.cm", "Rfam.clanin", "Rfam_metadata.tsv"):
        shutil.copy(os.path.join(RFAM_DIR, name), db / name)
    subprocess.run(
        ["cmpress", "-F", str(db / "Rfam.cm")], check=True, capture_output=True
    )
    return db


@pytest.mark.skipif(
    shutil.which("cmscan") is None or shutil.which("cmpress") is None,
    reason="Infernal not installed",
)
class TestCmscanIntegration:
    """End-to-end against real Infernal, using the 5-model test subset."""

    def _cmscan(self, db, fasta, out):
        subprocess.run(
            [
                "cmscan",
                "--rfam",
                "--cut_ga",
                "--nohmmonly",
                "--noali",
                "--fmt",
                "2",
                "--clanin",
                str(db / "Rfam.clanin"),
                "--tblout",
                str(out),
                str(db / "Rfam.cm"),
                fasta,
            ],
            check=True,
            capture_output=True,
        )
        return str(out)

    def test_finds_the_three_srnas_in_NC_004617(self, pressed_db, tmp_path):
        tblout = self._cmscan(pressed_db, NC_004617, tmp_path / "out.tblout")
        df = parse_cmscan_tblout(tblout, metadata=load_rfam_metadata(str(pressed_db)))

        assert sorted(df["rfam_id"].to_list()) == ["SprD", "SprX", "rli28"]
        assert set(df["type"].to_list()) == {"Gene; sRNA;"}

    def test_rfam_does_not_find_the_aragorn_tmrna(self, pressed_db, tmp_path):
        """Regression lock for a measured result.

        NC_051700 is pharokka's tmRNA test case and ARAGORN calls a tmRNA in
        it, but Rfam's RF00023 does not hit at all.  Rfam must therefore stay
        additive - it can never replace ARAGORN or tRNAscan-SE.
        """
        tblout = self._cmscan(pressed_db, NC_051700, tmp_path / "out.tblout")
        df = parse_cmscan_tblout(tblout, keep_trna=True)

        assert "RF00023" not in df["rfam_acc"].to_list()
