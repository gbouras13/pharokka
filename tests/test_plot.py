"""Tests for the circular genome plot, in particular ncRNA drawing.

Two genomes drove the design and both are in the test data:

``NC_004617.1``
    three Rfam sRNAs, every one of them wholly inside a CDS - the case that
    decided ncRNAs get a ring of their own rather than being painted into the
    CDS ring the way tRNA/tmRNA/CRISPR are.
``NC_028808.2``
    four vrRNA guide arrays inside 1.6 kb, all carrying the identical product
    "vrRNA guide array" - the case that decided the labels have to be thinned.

The ``golden/rfam`` fixture is NC_004617.1 run through pharokka, so the Rfam
half of this is exercised against a real pharokka GFF and GenBank rather than
a hand-written one.  vrRNA arrays are appended to a copy of that GFF, because
none of the golden genomes carries one.
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import pytest
from Bio.SeqFeature import FeatureLocation, SeqFeature
from pycirclize import Circos
from pycirclize.parser import Genbank

matplotlib.use("Agg")

from pharokka.plot import (  # noqa: E402
    NCRNA_COL,
    NCRNA_TRACK_RANGE,
    PLOT_FEATURE_TYPES,
    VRRNA_COL,
    VRRNA_PRODUCT,
    add_ncrna_track,
    create_plot,
    create_single_plot,
    is_vrrna_feature,
    ncrna_labels,
    ncrna_legend_handles,
    split_ncrna_features,
)

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")
GOLDEN_RFAM = os.path.join(TEST_DATA, "golden", "rfam")
GOLDEN_CRISPR = os.path.join(TEST_DATA, "golden", "crispr")

# the four NC_028808.2 arrays, as pharokka writes them
VRRNA_GFF_ROWS = [
    "NC_004617.1\tab initio prediction:pharokka:1.11.0\tncRNA\t27621\t27670\t26.0\t-\t.\t"
    "ID=PHARTEST_vrRNA_1;locus_tag=PHARTEST_vrRNA_1;product=vrRNA guide array;"
    "ncRNA_class=other;note=VIPR vrRNA guide array (strict call) | 10 GGYNN units\n",
    "NC_004617.1\tab initio prediction:pharokka:1.11.0\tncRNA\t27749\t27818\t50.7\t-\t.\t"
    "ID=PHARTEST_vrRNA_2;locus_tag=PHARTEST_vrRNA_2;product=vrRNA guide array;"
    "ncRNA_class=other;note=VIPR vrRNA guide array (strict call) | 14 GGYNN units\n",
    "NC_004617.1\tab initio prediction:pharokka:1.11.0\tncRNA\t27892\t27966\t86.0\t-\t.\t"
    "ID=PHARTEST_vrRNA_3;locus_tag=PHARTEST_vrRNA_3;product=vrRNA guide array;"
    "ncRNA_class=other;note=VIPR vrRNA guide array (strict call) | 15 GGYNN units\n",
    "NC_004617.1\tab initio prediction:pharokka:1.11.0\tncRNA\t29172\t29226\t121.0\t-\t.\t"
    "ID=PHARTEST_vrRNA_4;locus_tag=PHARTEST_vrRNA_4;product=vrRNA guide array;"
    "ncRNA_class=other;note=VIPR vrRNA guide array (strict call) | 11 GGYNN units\n",
]


def ncrna_feature(start, end, product, dbxref=None, strand=1):
    """A bare ncRNA SeqFeature, shaped like the ones both parsers hand back."""
    qualifiers = {"product": [product]}
    if dbxref is not None:
        qualifiers["Dbxref"] = [dbxref]
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="ncRNA",
        qualifiers=qualifiers,
    )


def rfam_feature(start, end, product="Cobalamin riboswitch", acc="RF00174", **kwargs):
    return ncrna_feature(start, end, product, dbxref=f"RFAM:{acc}", **kwargs)


def vrrna_feature(start, end, **kwargs):
    return ncrna_feature(start, end, VRRNA_PRODUCT, **kwargs)


def write_gff(tmp_path, extra_rows=(), drop_ncrna=False):
    """A copy of the golden rfam GFF, optionally edited."""
    path = tmp_path / "pharokka.gff"
    with open(os.path.join(GOLDEN_RFAM, "pharokka.gff")) as fin:
        rows = fin.readlines()
    if drop_ncrna:
        rows = [row for row in rows if "\tncRNA\t" not in row]
    with open(path, "w") as fout:
        fout.writelines(rows)
        fout.writelines(extra_rows)
    return str(path)


def plot_to(tmp_path, gff_file, gbk_file=None, **kwargs):
    """Run create_plot with pharokka's own defaults and return the SVG text."""
    png = tmp_path / "plot.png"
    svg = tmp_path / "plot.svg"
    options = dict(
        interval=5000,
        annotations=1,
        title_size=20,
        plot_title="test",
        truncate=20,
        dpi=100,
        label_size=8,
        label_hypotheticals=False,
        remove_other_features_labels=False,
        label_force_list=[],
    )
    options.update(kwargs)
    create_plot(
        gff_file=gff_file,
        gbk_file=gbk_file or os.path.join(GOLDEN_RFAM, "pharokka.gbk"),
        outfile=str(png),
        svg_plot_file=str(svg),
        **options,
    )
    assert png.exists() and svg.exists()
    return svg.read_text()


class TestIsVrrnaFeature:
    """vrRNA arrays and Rfam hits are both ncRNA features; the attributes
    are the only thing telling them apart."""

    def test_an_rfam_dbxref_is_not_a_vrrna(self):
        assert is_vrrna_feature(rfam_feature(100, 250)) is False

    def test_the_vrrna_product_is_a_vrrna(self):
        assert is_vrrna_feature(vrrna_feature(100, 150)) is True

    def test_an_ncrna_with_neither_marker_is_not_a_vrrna(self):
        assert is_vrrna_feature(ncrna_feature(100, 150, "some other RNA")) is False

    def test_the_rfam_dbxref_wins_over_the_product(self):
        # an Rfam family described as a vrRNA guide array is still an Rfam hit
        feature = rfam_feature(100, 150, product=VRRNA_PRODUCT)
        assert is_vrrna_feature(feature) is False

    def test_a_non_rfam_dbxref_does_not_make_it_an_rfam_hit(self):
        feature = ncrna_feature(100, 150, VRRNA_PRODUCT, dbxref="SO:0000655")
        assert is_vrrna_feature(feature) is True


class TestSplitNcrnaFeatures:
    def test_splits_into_rfam_and_vrrna(self):
        features = [
            rfam_feature(100, 250),
            vrrna_feature(500, 550),
            rfam_feature(900, 1050),
        ]
        rfam, vrrna = split_ncrna_features(features)
        assert len(rfam) == 2
        assert len(vrrna) == 1
        assert int(vrrna[0].location.start) == 500

    def test_empty_input(self):
        assert split_ncrna_features([]) == ([], [])


class TestNcrnaLabels:
    def test_the_product_is_the_label(self):
        pos_list, labels = ncrna_labels([rfam_feature(100, 250)], truncate=40)
        assert labels == ["Cobalamin riboswitch"]
        assert pos_list == [175.0]

    def test_long_labels_are_truncated(self):
        feature = rfam_feature(100, 250, product="Small pathogenicity island RNA X")
        _, labels = ncrna_labels([feature], truncate=20)
        assert labels == ["Small pathogenicity ..."]

    def test_a_feature_without_a_product_is_skipped(self):
        feature = SeqFeature(
            FeatureLocation(100, 250, strand=1), type="ncRNA", qualifiers={}
        )
        assert ncrna_labels([feature], truncate=20) == ([], [])

    def test_empty_input(self):
        assert ncrna_labels([], truncate=20) == ([], [])

    def test_the_vrrna_cluster_is_thinned(self):
        """The four NC_028808.2 arrays sit inside 1.6 kb and all say
        "vrRNA guide array"; unthinned that is four labels on top of each
        other.  The two within 500 bp of the one before them go."""
        features = [
            vrrna_feature(27621, 27670),
            vrrna_feature(27749, 27818),
            vrrna_feature(27892, 27966),
            vrrna_feature(29172, 29226),
        ]
        pos_list, labels = ncrna_labels(features, truncate=20)
        assert labels == [VRRNA_PRODUCT, VRRNA_PRODUCT]
        assert pos_list == [pytest.approx(27645.5), pytest.approx(29199.0)]

    def test_well_spaced_features_are_all_kept(self):
        """The three NC_004617.1 sRNAs are kilobases apart and all survive."""
        features = [
            rfam_feature(36836, 37015, product="Listeria sRNA rli28"),
            rfam_feature(39255, 39405, product="Small pathogenicity island RNA X"),
            rfam_feature(40870, 41011, product="Small pathogenicity island RNA D"),
        ]
        pos_list, labels = ncrna_labels(features, truncate=40)
        assert len(pos_list) == 3
        assert labels[0] == "Listeria sRNA rli28"

    def test_both_kinds_are_thinned_together(self):
        """An array next to a riboswitch must not produce two labels on top
        of each other, so the two classes are thinned as one list."""
        features = [rfam_feature(1000, 1150), vrrna_feature(1200, 1250)]
        _, labels = ncrna_labels(features, truncate=40)
        assert labels == ["Cobalamin riboswitch"]

    def test_input_is_sorted_before_thinning(self):
        """The GFF holds the Rfam and vrRNA blocks separately, so the merged
        list arrives out of order."""
        features = [vrrna_feature(29172, 29226), vrrna_feature(27621, 27670)]
        pos_list, _ = ncrna_labels(features, truncate=20)
        assert pos_list == sorted(pos_list)
        assert len(pos_list) == 2

    def test_truncate_is_accepted_as_a_string(self):
        # create_plot has not yet cast its CLI argument at this point
        _, labels = ncrna_labels([rfam_feature(100, 250)], truncate="10")
        assert labels == ["Cobalamin ..."]


class TestNcrnaLegendHandles:
    def test_no_handles_without_ncrnas(self):
        assert ncrna_legend_handles([], []) == []

    def test_one_handle_per_kind_present(self):
        handles = ncrna_legend_handles([rfam_feature(1, 2)], [])
        assert [h.get_label() for h in handles] == ["ncRNA"]

        handles = ncrna_legend_handles([], [vrrna_feature(1, 2)])
        assert [h.get_label() for h in handles] == [VRRNA_PRODUCT]

        handles = ncrna_legend_handles([rfam_feature(1, 2)], [vrrna_feature(3, 4)])
        assert [h.get_label() for h in handles] == ["ncRNA", VRRNA_PRODUCT]


class TestAddNcrnaTrack:
    def sector(self):
        return Circos(sectors={"test": 50000}).get_sector("test")

    def test_no_track_is_added_without_ncrnas(self):
        """Genomes without ncRNAs must plot exactly as they did before."""
        sector = self.sector()
        assert add_ncrna_track(sector, [], []) is None
        assert sector.tracks == []

    def test_a_track_is_added_for_ncrnas(self):
        sector = self.sector()
        track = add_ncrna_track(sector, [rfam_feature(100, 250)], [])
        assert track is not None
        assert track.r_lim == NCRNA_TRACK_RANGE
        assert sector.tracks == [track]

    def test_a_track_is_added_for_vrrnas_alone(self):
        sector = self.sector()
        assert add_ncrna_track(sector, [], [vrrna_feature(100, 150)]) is not None

    def test_the_track_sits_inside_the_cds_track(self):
        # the CDS track is (70, 80) and GC content starts at 60
        assert NCRNA_TRACK_RANGE[1] <= 70
        assert NCRNA_TRACK_RANGE[0] >= 60


class TestPlotFeatureTypes:
    """pycirclize's GenBank parser defaults to CDS only, which is why the
    multiplotter used to silently drop every tRNA, tmRNA, CRISPR and ncRNA."""

    def test_every_type_create_single_plot_draws_is_requested(self):
        for feature_type in ("CDS", "tRNA", "tmRNA", "repeat_region", "ncRNA"):
            assert feature_type in PLOT_FEATURE_TYPES

    def test_ncrnas_reach_the_multiplotter(self):
        gbk = Genbank(os.path.join(GOLDEN_RFAM, "pharokka.gbk"))
        features = gbk.get_seqid2features(feature_type=PLOT_FEATURE_TYPES)
        ncrnas = [f for fs in features.values() for f in fs if f.type == "ncRNA"]
        assert len(ncrnas) == 3
        assert gbk.get_seqid2features() == {
            seqid: [f for f in fs if f.type == "CDS"] for seqid, fs in features.items()
        }

    def test_crisprs_reach_the_multiplotter(self):
        gbk = Genbank(os.path.join(GOLDEN_CRISPR, "pharokka.gbk"))
        features = gbk.get_seqid2features(feature_type=PLOT_FEATURE_TYPES)
        crisprs = [
            f for fs in features.values() for f in fs if f.type == "repeat_region"
        ]
        assert len(crisprs) == 2


class TestCreatePlot:
    """End to end over the golden NC_004617.1 GFF and GenBank."""

    def test_rfam_ncrnas_are_drawn(self, tmp_path):
        svg = plot_to(tmp_path, write_gff(tmp_path))
        assert NCRNA_COL.lower()[1:] in svg
        assert VRRNA_COL.lower()[1:] not in svg

    def test_vrrna_arrays_are_drawn_and_coloured_apart(self, tmp_path):
        svg = plot_to(tmp_path, write_gff(tmp_path, extra_rows=VRRNA_GFF_ROWS))
        assert NCRNA_COL.lower()[1:] in svg
        assert VRRNA_COL.lower()[1:] in svg

    def test_vrrna_arrays_alone(self, tmp_path):
        gff = write_gff(tmp_path, extra_rows=VRRNA_GFF_ROWS, drop_ncrna=True)
        svg = plot_to(tmp_path, gff)
        assert VRRNA_COL.lower()[1:] in svg
        assert NCRNA_COL.lower()[1:] not in svg

    def test_a_genome_without_ncrnas_is_unchanged(self, tmp_path):
        svg = plot_to(tmp_path, write_gff(tmp_path, drop_ncrna=True))
        assert NCRNA_COL.lower()[1:] not in svg
        assert VRRNA_COL.lower()[1:] not in svg

    def test_ncrnas_are_still_drawn_without_their_labels(self, tmp_path):
        """--remove_other_features_labels drops the labels, not the features."""
        svg = plot_to(tmp_path, write_gff(tmp_path), remove_other_features_labels=True)
        assert NCRNA_COL.lower()[1:] in svg

    def test_more_tmrnas_than_crisprs_does_not_raise(self, tmp_path):
        """The CRISPR label thinning used to loop over the *tmRNA* count and
        index into the CRISPR lists, so a genome with more tmRNAs than
        CRISPRs raised IndexError."""
        rows = [
            "NC_004617.1\tAragorn\ttmRNA\t1000\t1300\t.\t+\t.\tID=PHARTEST_tmRNA_1\n",
            "NC_004617.1\tAragorn\ttmRNA\t5000\t5300\t.\t+\t.\tID=PHARTEST_tmRNA_2\n",
            "NC_004617.1\tAragorn\ttmRNA\t9000\t9300\t.\t+\t.\tID=PHARTEST_tmRNA_3\n",
            "NC_004617.1\tminced\trepeat_region\t20000\t20100\t.\t+\t.\tID=PHARTEST_CRISPR_1\n",
        ]
        plot_to(tmp_path, write_gff(tmp_path, extra_rows=rows))


class TestCreateSinglePlot:
    """The multiplotter reads GenBank, so its feature filter is its own."""

    def single_plot(self, tmp_path, gbk_path, **kwargs):
        gbk = Genbank(gbk_path)
        seq_dict = gbk.get_seqid2seq()
        contig_id, contig_sequence = next(iter(seq_dict.items()))
        options = dict(
            interval=5000,
            annotations=1.0,
            title_size=20,
            plot_title="test",
            truncate=20,
            dpi=100,
            label_size=8,
            label_hypotheticals=False,
            remove_other_features_labels=False,
            label_force_list=[],
        )
        options.update(kwargs)
        create_single_plot(
            contig_id=contig_id,
            contig_sequence=contig_sequence,
            contig_count=len(seq_dict),
            gb_size_dict=gbk.get_seqid2size(),
            gb_feature_dict=gbk.get_seqid2features(feature_type=PLOT_FEATURE_TYPES),
            gbk=gbk,
            output=str(tmp_path),
            **options,
        )
        svg = tmp_path / f"{contig_id}.svg"
        assert (tmp_path / f"{contig_id}.png").exists()
        assert svg.exists()
        return svg.read_text()

    def test_rfam_ncrnas_are_drawn(self, tmp_path):
        svg = self.single_plot(tmp_path, os.path.join(GOLDEN_RFAM, "pharokka.gbk"))
        assert NCRNA_COL.lower()[1:] in svg
        assert VRRNA_COL.lower()[1:] not in svg

    def test_vrrna_arrays_are_coloured_apart(self, tmp_path):
        gbk_path = os.path.join(GOLDEN_RFAM, "pharokka.gbk")
        gbk = Genbank(gbk_path)
        feature_dict = gbk.get_seqid2features(feature_type=PLOT_FEATURE_TYPES)
        contig_id = next(iter(gbk.get_seqid2seq()))
        feature_dict[contig_id].append(vrrna_feature(27621, 27670, strand=-1))

        seq_dict = gbk.get_seqid2seq()
        create_single_plot(
            contig_id=contig_id,
            contig_sequence=seq_dict[contig_id],
            contig_count=1,
            gb_size_dict=gbk.get_seqid2size(),
            gb_feature_dict=feature_dict,
            gbk=gbk,
            interval=5000,
            annotations=1.0,
            title_size=20,
            plot_title="test",
            truncate=20,
            output=str(tmp_path),
            dpi=100,
            label_size=8,
            label_hypotheticals=False,
            remove_other_features_labels=False,
            label_force_list=[],
        )
        svg = (tmp_path / f"{contig_id}.svg").read_text()
        assert NCRNA_COL.lower()[1:] in svg
        assert VRRNA_COL.lower()[1:] in svg

    def test_a_genome_without_ncrnas_draws_no_ncrna_track(self, tmp_path):
        svg = self.single_plot(tmp_path, os.path.join(GOLDEN_CRISPR, "pharokka.gbk"))
        assert NCRNA_COL.lower()[1:] not in svg
        assert VRRNA_COL.lower()[1:] not in svg


@pytest.fixture(autouse=True)
def close_figures():
    """create_plot leaves its matplotlib figure open."""
    yield
    plt.close("all")


def test_golden_fixtures_are_present():
    """Guards the two fixture directories these tests read."""
    for path in (GOLDEN_RFAM, GOLDEN_CRISPR):
        assert os.path.isfile(os.path.join(path, "pharokka.gff"))
        assert os.path.isfile(os.path.join(path, "pharokka.gbk"))
