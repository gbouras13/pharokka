"""Tests for the circular genome plot.

The fixtures are pharokka's own golden outputs: ``golden/rfam`` is
NC_004617.1 and ``golden/crispr`` is the CRISPR case, so the parsing these
tests do is the parsing a real plot does.
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import pytest
from pycirclize.parser import Genbank

matplotlib.use("Agg")

from pharokka.plot import PLOT_FEATURE_TYPES, create_plot  # noqa: E402

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")
GOLDEN_RFAM = os.path.join(TEST_DATA, "golden", "rfam")
GOLDEN_CRISPR = os.path.join(TEST_DATA, "golden", "crispr")


def write_gff(tmp_path, extra_rows=()):
    """A copy of the golden rfam GFF with extra feature rows appended."""
    path = tmp_path / "pharokka.gff"
    with open(os.path.join(GOLDEN_RFAM, "pharokka.gff")) as fin:
        rows = fin.readlines()
    with open(path, "w") as fout:
        fout.writelines(rows)
        fout.writelines(extra_rows)
    return str(path)


def plot_to(tmp_path, gff_file, **kwargs):
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
        gbk_file=os.path.join(GOLDEN_RFAM, "pharokka.gbk"),
        outfile=str(png),
        svg_plot_file=str(svg),
        **options,
    )
    assert png.exists() and svg.exists()
    return svg.read_text()


class TestPlotFeatureTypes:
    """pycirclize's GenBank parser defaults to CDS only, which is why the
    multiplotter used to silently drop every tRNA, tmRNA and CRISPR."""

    def test_every_type_create_single_plot_draws_is_requested(self):
        for feature_type in ("CDS", "tRNA", "tmRNA", "repeat_region"):
            assert feature_type in PLOT_FEATURE_TYPES

    def test_the_pycirclize_default_really_is_cds_only(self):
        gbk = Genbank(os.path.join(GOLDEN_RFAM, "pharokka.gbk"))
        features = gbk.get_seqid2features(feature_type=PLOT_FEATURE_TYPES)
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


class TestCrisprLabelThinning:
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
