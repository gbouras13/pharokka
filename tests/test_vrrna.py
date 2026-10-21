"""Tests for de novo vrRNA guide array detection.

The three golden genomes are the ones the design was settled on:

``NC_028808.2``
    *Escherichia* phage SUSP1, the type locus of Lopez *et al.* 2026 - four
    arrays in 1.6 kb, two of them abutting, two tail families.
``NC_074658.1``
    *Pseudomonas* phage Medea1, a blind validation: one array, found by the
    caller before anyone had looked at this genome for VIPR.
``NC_073440.1`` and ``CP000623.1``
    high-GC negatives, 66.3% and 65.4% GC, which is where the naive caller
    produces nothing but false positives.  They fail at different stages and
    both stages matter: Albright never gets past the periodicity prefilter,
    phiE202 does and is then rejected by the shuffled null.
"""

import os

import polars as pl
import pytest

from pharokka.post_processing import Pharok
from pharokka.vrrna import (
    CLUSTER_WINDOW_BP,
    DEFAULT_TAIL,
    add_locus_tags,
    dinuc_shuffle,
    empirical_p,
    empty_vrrna_df,
    fixed_matches,
    null_max_units,
    periodicity_windows,
    resolve,
    revcomp,
    scan,
    scan_contig,
    scan_fasta,
    tail_identity,
    unambiguous_nn,
    unit_ok,
    vrrna_note,
    write_vrrna_tsv,
)

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")
VRRNA_DIR = os.path.join(TEST_DATA, "vrrna")

SUSP1 = os.path.join(VRRNA_DIR, "NC_028808.2.fasta")
MEDEA1 = os.path.join(VRRNA_DIR, "NC_074658.1.fasta")
ALBRIGHT = os.path.join(VRRNA_DIR, "NC_073440.1.fasta")
PHIE202 = os.path.join(VRRNA_DIR, "CP000623.1.fasta")


def _seq(path):
    return "".join(
        line.strip() for line in open(path) if not line.startswith(">")
    ).upper()


@pytest.fixture(scope="module")
def susp1():
    return _seq(SUSP1)


@pytest.fixture(scope="module")
def medea1():
    return _seq(MEDEA1)


@pytest.fixture(scope="module")
def albright():
    return _seq(ALBRIGHT)


@pytest.fixture(scope="module")
def phie202():
    return _seq(PHIE202)


# The shuffled null re-scans a whole genome 100 times, so it is seconds per
# genome rather than milliseconds.  Everything below that needs a scored call
# set shares one, module-scoped - the caller is deterministic, so a cached
# result is the same result.


@pytest.fixture(scope="module")
def susp1_calls(susp1):
    return scan_contig(susp1)


@pytest.fixture(scope="module")
def medea1_calls(medea1):
    return scan_contig(medea1)


@pytest.fixture(scope="module")
def susp1_df():
    return scan_fasta(SUSP1)


@pytest.fixture(scope="module")
def medea1_df():
    return scan_fasta(MEDEA1)


class TestUnits:
    @pytest.mark.parametrize("unit", ["GGCAT", "GGTCT", "GGCGG", "GGTTT"])
    def test_conforming_units(self, unit):
        assert unit_ok(unit)
        assert fixed_matches(unit) == 3

    @pytest.mark.parametrize(
        "unit,matches",
        [("GGAAT", 2), ("GCCAT", 2), ("AGCAT", 2), ("AACAT", 1), ("AAAAT", 0)],
    )
    def test_non_conforming_units(self, unit, matches):
        assert not unit_ok(unit)
        assert fixed_matches(unit) == matches

    def test_y_is_c_or_t_only(self):
        assert unit_ok("GGCAT") and unit_ok("GGTAT")
        assert not unit_ok("GGAAT") and not unit_ok("GGGAT")

    def test_wrong_length_never_conforms(self):
        assert not unit_ok("GGCA")
        assert fixed_matches("GGCA") == 0

    def test_ambiguous_guide_positions_rejected(self):
        assert unambiguous_nn("GGCAT")
        assert not unambiguous_nn("GGCNT")

    def test_revcomp_round_trips(self):
        assert revcomp(revcomp("GGCATGGTCTAAN")) == "GGCATGGTCTAAN"


class TestPeriodicityPrefilter:
    def test_finds_a_tandem_array(self):
        flank = "ACTGAACCTTGACATGCATCAGTTACGGAATCCGTAAGCTA" * 10
        array = "GGCATGGTCTGGTCTGGTGTGGTGAGGTATGGCCAGGTTTGGCGCGGTAA" * 2
        seq = flank + array + flank
        windows = periodicity_windows(seq)
        assert windows, "a 100 nt period-5 array must survive the prefilter"
        assert any(s <= len(flank) and e >= len(flank) + len(array) for s, e in windows)

    def test_short_sequence_is_passed_through_whole(self):
        assert periodicity_windows("ACGT" * 5) == [(0, 20)]

    def test_density_floor_is_what_short_arrays_fail(self):
        # 0.45 rather than 0.55: a 6-unit array is 30 nt and cannot fill the
        # 40 nt window on its own.  Raising the floor loses it.
        flank = "ACTGAACCTTGACATGCATCAGTTACGGAATCCGTAAGCTA" * 10
        seq = flank + "GGCATGGTCTGGTCTGGTGTGGTGAGGTAT" + flank
        assert periodicity_windows(seq, density=0.45) != []


class TestBlockConformity:
    def test_extension_bridges_an_early_bad_unit(self):
        # From RECALL_BENCHMARK.md 4.1.  This block is 0.810 conforming
        # overall but only 0.667 after two units, so a prefix-conformity rule
        # stops dead at unit 2 and loses the array; block conformity keeps it.
        units = ["GGCGC", "GCGGC", "GGCGC", "GGTCA", "TACCC", "GGTAG", "GGCCG"]
        block = "".join(units)
        assert sum(fixed_matches(u) for u in units) / (3 * len(units)) == pytest.approx(
            0.810, abs=0.001
        )
        # the register walk stops 5 nt short of the window end, so the block
        # needs a non-conforming tail for its own last unit to be in scope
        calls = scan(block + "AAAAAAAAAA", min_units=6, use_prefilter=False)
        assert [(c["start"], c["end"], c["units"]) for c in calls] == [(0, 35, 7)]

    def test_prefix_conformity_would_have_lost_it(self):
        # the same block truncated where a prefix rule would stop: two units
        # in, conformity is 0.667 and below the 0.80 floor
        units = ["GGCGC", "GCGGC"]
        assert sum(fixed_matches(u) for u in units) / (3 * len(units)) == pytest.approx(
            0.667, abs=0.001
        )


class TestResolve:
    def _call(self, start, end, units, strand="+", conformity=0.9):
        return {
            "start": start,
            "end": end,
            "units": units,
            "strand": strand,
            "conformity": conformity,
        }

    def test_abutting_arrays_both_survive(self):
        # Adjacent tandem arrays in different registers share exactly 1 bp.
        # An any-overlap rule deletes one of every such pair; the SUSP1 locus
        # has two.
        kept = resolve([self._call(0, 51, 10), self._call(50, 120, 14)])
        assert len(kept) == 2

    def test_real_overlap_suppresses_the_smaller_call(self):
        kept = resolve([self._call(0, 60, 10), self._call(40, 110, 14)])
        assert len(kept) == 1
        assert kept[0]["units"] == 14

    def test_opposite_strands_never_suppress_each_other(self):
        kept = resolve([self._call(0, 60, 10, "+"), self._call(0, 60, 14, "-")])
        assert len(kept) == 2

    def test_seed_calls_survive_unconditionally(self):
        seed = [self._call(0, 60, 6)]
        kept = resolve([self._call(10, 80, 14)], seed=seed)
        assert kept == seed


class TestTail:
    def test_identical_tail_scores_one(self):
        assert tail_identity(DEFAULT_TAIL, DEFAULT_TAIL) == 1.0

    def test_short_tail_scores_zero(self):
        assert tail_identity("GGGTAGTTAG", DEFAULT_TAIL) == 0.0

    def test_unrelated_tail_scores_low(self):
        assert tail_identity("A" * 60, DEFAULT_TAIL) < 0.5


class TestNullModels:
    def test_mono_shuffle_preserves_base_composition(self, albright):
        sub = albright[:5000]
        import random

        from pharokka.vrrna import _mono_shuffler

        shuffled = next(_mono_shuffler(sub, random.Random(1)))
        assert sorted(shuffled) == sorted(sub)

    def test_dinuc_shuffle_preserves_dinucleotide_counts(self, albright):
        import collections
        import random

        sub = albright[:5000]
        shuffled = dinuc_shuffle(sub, random.Random(1))

        def dinucs(s):
            return collections.Counter(s[i : i + 2] for i in range(len(s) - 1))

        assert sorted(shuffled) == sorted(sub)
        assert dinucs(shuffled) == dinucs(sub)
        assert shuffled[0] == sub[0] and shuffled[-1] == sub[-1]
        assert shuffled != sub

    def test_null_is_deterministic_for_a_seed(self, albright):
        sub = albright[:8000]
        assert null_max_units(sub, 5, 6, seed=7) == null_max_units(sub, 5, 6, seed=7)

    def test_unknown_null_model_is_rejected(self, albright):
        with pytest.raises(ValueError, match="unknown null model"):
            null_max_units(albright[:1000], 1, 6, model="trinuc")

    @pytest.mark.parametrize(
        "units,nulls,expected",
        [
            (10, [0] * 99, 0.01),
            (6, [6] * 99, 1.0),
            (7, [6] * 49 + [7] * 50, 0.51),
        ],
    )
    def test_empirical_p(self, units, nulls, expected):
        assert empirical_p(units, nulls) == pytest.approx(expected, abs=0.005)

    def test_empirical_p_without_a_null_is_nan(self):
        assert empirical_p(10, []) != empirical_p(10, [])  # NaN


class TestGoldenSUSP1:
    """NC_028808.2 - the VIPR type locus.

    Note the first array.  Lopez *et al.* publish it as 27620-27650 (6 units);
    the block-conformity extension legitimately reads it as 27620-27670 (10
    units), because GGTCC conforms and the ten-unit block still clears the
    0.80 floor.  The paper's caller evidently applies a stricter test at the
    5' end, but the leading-window guard that reproduces it was measured and
    costs twelve arrays of boundary accuracy elsewhere to gain this one.  So
    the expectation here is the longer call, deliberately, and the caller is
    not bent to fit the paper.  RECALL_BENCHMARK.md 8.1 has the numbers.
    """

    EXPECTED = [
        (27620, 27670, "-", 10),
        (27748, 27818, "-", 14),
        (27891, 27966, "-", 15),
        (29171, 29226, "-", 11),
    ]

    def test_four_arrays_at_the_published_coordinates(self, susp1_calls):
        assert [
            (c["fstart"], c["fend"], c["strand"], c["units"]) for c in susp1_calls
        ] == self.EXPECTED

    def test_every_array_clears_the_null(self, susp1_calls):
        calls = susp1_calls
        assert all(c["pvalue"] <= 0.05 for c in calls)
        # 100 shuffles, so 1/101 is the floor and all four sit on it
        assert all(c["pvalue"] == pytest.approx(0.0099, abs=1e-4) for c in calls)

    def test_all_four_are_strict_calls(self, susp1_calls):
        assert {c["tier"] for c in susp1_calls} == {"strict"}

    def test_the_abutting_pair_both_survive(self, susp1_calls):
        # 27748-27818 and 27891-27966 are the pair that an any-overlap
        # resolution rule would collapse into one.
        spans = {(c["fstart"], c["fend"]) for c in susp1_calls}
        assert (27748, 27818) in spans and (27891, 27966) in spans

    def test_two_tail_families(self, susp1_calls):
        calls = susp1_calls
        exact = [c for c in calls if c["tail_id"] == 1.0]
        assert len(exact) == 2, "two arrays carry the type II consensus tail"
        assert all(c["tail_id"] < 0.7 for c in calls if c["tail_id"] != 1.0)

    def test_tail_agreement_is_reported_within_the_genome(self, susp1_calls):
        calls = susp1_calls
        by_span = {(c["fstart"], c["fend"]): c for c in calls}
        # the two type II arrays agree with each other far better than either
        # agrees with the other family
        assert by_span[(27748, 27818)]["tail_agreement"] > 0.8
        assert by_span[(27620, 27670)]["tail_agreement"] < 0.7

    def test_all_four_arrays_form_one_cluster(self, susp1_calls):
        calls = susp1_calls
        assert {c["cluster_id"] for c in calls} == {1}
        assert {c["cluster_size"] for c in calls} == {4}
        assert calls[-1]["fend"] - calls[0]["fstart"] < CLUSTER_WINDOW_BP


class TestGoldenMedea1:
    """NC_074658.1 - one array plus the full VIPR locus, found blind."""

    def test_one_array(self, medea1_calls):
        assert [
            (c["fstart"], c["fend"], c["strand"], c["units"]) for c in medea1_calls
        ] == [(8025, 8085, "-", 12)]

    def test_it_clears_the_null(self, medea1_calls):
        assert medea1_calls[0]["pvalue"] == pytest.approx(0.0099, abs=1e-4)

    def test_a_lone_array_has_no_tail_agreement(self, medea1_calls):
        call = medea1_calls[0]
        assert call["tail_agreement"] is None
        assert call["cluster_size"] == 1

    def test_the_relaxed_tier_does_not_invent_neighbours(self, medea1_calls):
        # Dropping the unit floor to 4 inside the locus turns up three short
        # runs at 58% GC.  Scored against a null drawn from the same
        # neighbourhood they come out at p 0.29-0.53 and are correctly
        # dropped; without that local null they would all be reported.
        assert len(medea1_calls) == 1


class TestHighGCNegatives:
    def test_albright_is_rejected_by_the_prefilter(self, albright):
        # 66.3% GC.  Nothing survives the periodicity prefilter at the 0.45
        # density floor, so the null is never even consulted.  Switch the
        # prefilter off and a 7-unit array appears - which is the whole
        # reason the prefilter lowers the false-call rate rather than just
        # making the scan faster.
        assert scan(albright) == []
        assert len(scan(albright, use_prefilter=False)) == 1
        assert scan_contig(albright) == []

    def test_phie202_is_rejected_by_the_null(self, phie202):
        # 65.4% GC.  Here a 7-unit array does clear the prefilter, so this is
        # the case that tests the null itself rather than the prefilter.
        raw = scan(phie202)
        assert len(raw) == 1 and raw[0]["units"] == 7
        assert scan_contig(phie202) == []

    def test_relaxing_max_p_lets_the_phie202_call_back_through(self, phie202):
        calls = scan_contig(phie202, max_p=1.0)
        strict = [c for c in calls if c["tier"] == "strict"]
        assert len(strict) == 1
        assert strict[0]["pvalue"] > 0.5
        # with the gate wide open the relaxed tier fires too, which is what
        # max_p = 1 asks for
        assert len(calls) > 1

    def test_max_p_is_the_knob_not_the_unit_count(self, phie202):
        # The phiE202 call has 7 units, one more than the default floor, so no
        # setting of min_units separates it from a real 7-unit array.  Only
        # the null does.
        assert scan(phie202, min_units=7) != []
        assert scan_contig(phie202, max_p=0.05) == []


class TestTwoTier:
    def test_relaxed_calls_can_never_change_a_strict_one(self, susp1, susp1_calls):
        strict = scan(susp1)
        both = susp1_calls
        strict_spans = {(c["fstart"], c["fend"]) for c in strict}
        assert strict_spans <= {(c["fstart"], c["fend"]) for c in both}

    def test_a_protein_hit_opens_a_relaxed_window_on_its_own(self, albright):
        # The hook for the Vipr HMM: a protein hit seeds the relaxed tier even
        # where no array cleared the genome-wide null.  Albright has a 7-unit
        # array at 17996-18031 that the prefilter hides; a protein hit beside
        # it brings the region back into scope.  Whether anything is then
        # reported is still the local null's decision.
        without = scan_contig(albright)
        with_hit = scan_contig(albright, protein_hits=[(17000, 17900)], max_p=1.0)
        assert without == []
        assert any(c["tier"] == "relaxed" for c in with_hit)

    def test_relaxed_tier_is_off_when_its_floor_is_not_lower(self, susp1):
        calls = scan_contig(susp1, relaxed_min_units=6)
        assert {c["tier"] for c in calls} == {"strict"}


class TestDataFrame:
    def test_empty_schema_is_stable(self):
        df = empty_vrrna_df()
        assert df.height == 0
        assert df.columns == [
            "contig",
            "locus_tag",
            "start",
            "stop",
            "strand",
            "units",
            "conformity",
            "score",
            "pvalue",
            "tier",
            "cluster_id",
            "cluster_size",
            "tail_identity",
            "tail_agreement",
            "gc",
            "guides",
            "tail",
        ]

    def test_coordinates_become_one_based_inclusive(self, susp1_df):
        df = susp1_df
        # the caller works 0-based half-open, like the published array table;
        # everything downstream of here is 1-based inclusive
        assert df["start"].to_list() == [27621, 27749, 27892, 29172]
        assert df["stop"].to_list() == [27670, 27818, 27966, 29226]
        assert (df["stop"] - df["start"] + 1).to_list() == [50, 70, 75, 55]

    def test_span_is_five_nucleotides_per_unit(self, susp1_df):
        df = susp1_df
        assert (df["stop"] - df["start"] + 1).to_list() == [
            5 * u for u in df["units"].to_list()
        ]

    def test_contig_and_gc_are_recorded(self, medea1_df):
        df = medea1_df
        assert df["contig"].to_list() == ["NC_074658.1"]
        assert df["gc"][0] == pytest.approx(0.581, abs=0.001)

    def test_empty_result_still_has_the_schema(self):
        df = scan_fasta(ALBRIGHT)
        assert df.height == 0
        assert df.schema == empty_vrrna_df().schema

    def test_skipping_the_null_reports_every_raw_call(self):
        # n_shuffles=0 is the "show me everything" escape hatch - no null, so
        # no gate on either tier.  The CLI refuses it (--vrrna_shuffles must
        # be >= 1); it exists for people calling the module directly.
        df = scan_fasta(PHIE202, n_shuffles=0)
        assert df.filter(pl.col("tier") == "strict").height == 1
        assert df["pvalue"].is_nan().all()

    def test_dinuc_null_is_available_as_an_option(self):
        df = scan_fasta(MEDEA1, null_model="dinuc")
        assert df.height == 1
        assert df["pvalue"][0] <= 0.05


class TestLocusTags:
    def _df(self, contigs):
        return pl.DataFrame({"contig": contigs, "locus_tag": [""] * len(contigs)})

    def test_single_contig_numbering(self):
        tagged = add_locus_tags(self._df(["c1", "c1", "c1"]), "ABCDEFGH", 1)
        assert tagged["locus_tag"].to_list() == [
            "ABCDEFGH_vrRNA_1",
            "ABCDEFGH_vrRNA_2",
            "ABCDEFGH_vrRNA_3",
        ]

    def test_multi_contig_numbering_restarts_per_contig(self):
        tagged = add_locus_tags(self._df(["c1", "c1", "c2"]), "ABCDEFGH", 2)
        assert tagged["locus_tag"].to_list() == [
            "c1_vrRNA_1",
            "c1_vrRNA_2",
            "c2_vrRNA_1",
        ]

    def test_empty_frame_is_returned_unchanged(self):
        df = empty_vrrna_df()
        assert add_locus_tags(df, "ABCDEFGH", 1).height == 0


class TestNote:
    @pytest.fixture
    def row(self, susp1_df):
        return susp1_df.to_dicts()[1]

    def test_note_has_no_gff_structural_characters(self, row):
        note = vrrna_note(row)
        assert ";" not in note and "=" not in note and "," not in note

    def test_note_reports_the_p_value_and_the_unit_count(self, row):
        note = vrrna_note(row)
        assert "14 GGYNN units" in note
        assert "empirical p 0.0099" in note

    def test_note_describes_the_cluster_without_multiplying_it(self, row):
        note = vrrna_note(row)
        assert "4 arrays within 5000 bp" in note
        assert "correlated evidence" in note

    def test_note_reports_both_tail_measures(self, row):
        note = vrrna_note(row)
        assert "tail identity to type II consensus 1.000" in note
        assert "another array on this contig" in note

    def test_lone_array_note_omits_tail_agreement(self, medea1_df):
        note = vrrna_note(medea1_df.to_dicts()[0])
        assert "another array on this contig" not in note
        assert "cluster" not in note


class TestTsv:
    def test_written_file_round_trips(self, tmp_path, susp1_df):
        df = add_locus_tags(susp1_df, "PHARTEST", 1)
        write_vrrna_tsv(df, str(tmp_path), "pharokka")
        out = tmp_path / "pharokka_vrrna.tsv"
        assert out.exists()
        reread = pl.read_csv(out, separator="\t")
        assert reread.height == 4
        assert reread["locus_tag"].to_list() == [
            f"PHARTEST_vrRNA_{i}" for i in range(1, 5)
        ]

    def test_empty_frame_writes_a_header_only_file(self, tmp_path):
        write_vrrna_tsv(empty_vrrna_df(), str(tmp_path), "pharokka")
        lines = (tmp_path / "pharokka_vrrna.tsv").read_text().splitlines()
        assert len(lines) == 1
        assert lines[0].startswith("contig\tlocus_tag\tstart\tstop\tstrand")


class TestGffRows:
    def _pharok(self, susp1_df):
        pharok = Pharok()
        pharok.vrrna_df = susp1_df.clone()
        pharok.vrrna_flag = True
        pharok.locustag = "PHARTEST"
        pharok.length_df = pl.DataFrame({"contig": ["NC_028808.2"], "length": [90743]})
        return pharok

    def test_gff_columns_and_feature_type(self, susp1_df):
        gff = self._pharok(susp1_df)._build_vrrna_gff_df()
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
        assert gff["Region"].to_list() == ["ncRNA"] * 4
        assert gff["strand"].to_list() == ["-"] * 4
        assert all(
            m.startswith("ab initio prediction:pharokka:") for m in gff["Method"]
        )

    def test_attributes_are_well_formed(self, susp1_df):
        gff = self._pharok(susp1_df)._build_vrrna_gff_df()
        for attrs in gff["attributes"]:
            keys = [pair.split("=", 1)[0] for pair in attrs.split(";")]
            assert keys == ["ID", "locus_tag", "product", "ncRNA_class", "note"]
            assert "product=vrRNA guide array" in attrs
            assert "ncRNA_class=other" in attrs

    def test_locus_tags_are_written_back_to_the_frame(self, susp1_df):
        pharok = self._pharok(susp1_df)
        pharok._build_vrrna_gff_df()
        assert pharok.vrrna_df["locus_tag"].to_list() == [
            f"PHARTEST_vrRNA_{i}" for i in range(1, 5)
        ]

    def test_coordinates_match_the_dataframe(self, susp1_df):
        pharok = self._pharok(susp1_df)
        gff = pharok._build_vrrna_gff_df()
        assert gff["start"].to_list() == pharok.vrrna_df["start"].to_list()
        assert gff["stop"].to_list() == pharok.vrrna_df["stop"].to_list()
