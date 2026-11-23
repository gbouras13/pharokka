"""
Unit tests for pharokka

Usage: pytest

"""

# import
import subprocess as sp
import sys
import unittest
from pathlib import Path
from unittest import mock

import pytest
from loguru import logger

from pharokka.input_commands import (
    check_dependencies,
    instantiate_dirs,
    parse_mmseqs_version,
    validate_fasta,
    validate_gene_predictor,
    validate_meta,
    validate_strand,
    validate_terminase,
    validate_terminase_start,
    validate_threads,
)
from pharokka.util import remove_directory

# test data
test_data = Path("tests/test_data")
functions_data = Path(f"{test_data}/functions_files")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{overall_data}/Meta_example")
standard_data = Path(f"{overall_data}/Standard_examples")
logger.add(lambda _: sys.exit(1), level="ERROR")


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class TestInstantiateDir(unittest.TestCase):
    """Tests of instantiate dir"""

    def test_instantiate_dirs(self):
        fake_out_dir: Path = f"{test_data}/fake_out"
        instantiate_dirs(fake_out_dir, meta=False, force=False)
        remove_directory(fake_out_dir)

    def test_instantiate_overwrite(self):
        fake_out_dir_for_test: Path = f"{test_data}/fake_out_for_test"
        instantiate_dirs(fake_out_dir_for_test, meta=False, force=True)

    def test_instantiate_overwrite_fail(self):
        with self.assertRaises(SystemExit):
            fake_out_dir_for_test: Path = f"{test_data}/fake_out_for_test"
            instantiate_dirs(fake_out_dir_for_test, meta=False, force=False)

    def test_instantiate_dirs_fail_output_exists(self):
        with self.assertRaises(SystemExit):
            instantiate_dirs(test_data, meta=False, force=False)


class TestValidateFasta(unittest.TestCase):
    """Tests of validate_fasta"""

    def test_validate_fasta(self):
        fasta: Path = f"{functions_data}/good.fasta"
        validate_fasta(fasta)

    def test_validate_fasta_bad(self):
        with self.assertRaises(SystemExit):
            fasta: Path = f"{functions_data}/bad.fasta"
            validate_fasta(fasta)


class TestValidateGenePredictor(unittest.TestCase):
    """Tests of validate_gene_predictor"""

    def test_validate_gene_predictor(self):
        genbank_flag = False
        validate_gene_predictor("phanotate", genbank_flag)
        validate_gene_predictor("prodigal", genbank_flag)

    def test_validate_genbank_flag(self):
        genbank_flag = True
        validate_gene_predictor("genbank", genbank_flag)

    def test_gene_predictor_bad(self):
        genbank_flag = False
        with self.assertRaises(SystemExit):
            validate_gene_predictor("blah", genbank_flag)


class TestValidateMeta(unittest.TestCase):
    """Tests of validate_meta"""

    def test_validate_meta_good(self):
        fasta: Path = f"{meta_data}/fake_meta.fa"
        validate_meta(fasta, meta=True, split=False, genbank=False)

    def test_validate_meta_good_split(self):
        fasta: Path = f"{meta_data}/fake_meta.fa"
        validate_meta(fasta, meta=True, split=True, genbank=False)

    def test_validate_meta_bad(self):
        with self.assertRaises(SystemExit):
            fasta: Path = f"{standard_data}/SAOMS1.fasta"
            validate_meta(fasta, meta=True, split=True, genbank=False)

    def test_validate_no_meta_split(self):
        fasta: Path = f"{meta_data}/fake_meta.fa"
        validate_meta(fasta, meta=False, split=True, genbank=False)

    def test_validate_no_meta_no_split(self):
        fasta: Path = f"{meta_data}/fake_meta.fa"
        validate_meta(fasta, meta=False, split=False, genbank=False)


class TestValidateStrand(unittest.TestCase):
    """Tests of validate_strand"""

    def test_validate_strand(self):
        validate_strand("pos")
        validate_strand("neg")

    def test_validate_strand_bad(self):
        with self.assertRaises(SystemExit):
            validate_strand("blah")


class TestValidateTerminaseStart(unittest.TestCase):
    """Tests of validate_terminase_start"""

    def test_validate_terminase_start(self):
        validate_terminase_start(10)

    def test_validate_terminase_start_bad(self):
        with self.assertRaises(SystemExit):
            validate_terminase_start("blah")


class TestValidateTerminase(unittest.TestCase):
    """Tests of validate_terminase"""

    def test_validate_terminase(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        validate_terminase(fasta, "pos", 323)

    def test_validate_terminase_bad_multiple_contigs(self):
        with self.assertRaises(SystemExit):
            fasta: Path = f"{meta_data}/fake_meta.fa"
            validate_terminase(fasta, "pos", 323)

    def test_validate_terminase_bad_nothing(self):
        with self.assertRaises(SystemExit):
            fasta: Path = f"{meta_data}/fake_meta.fa"
            validate_terminase(fasta, "nothing", 323)

    def test_validate_terminase_bad_nothing2(self):
        with self.assertRaises(SystemExit):
            fasta: Path = f"{meta_data}/fake_meta.fa"
            validate_terminase(fasta, "pos", "nothing")


class TestValidateThreads(unittest.TestCase):
    """Tests of validate_threads"""

    def test_validate_threads(self):
        validate_threads(1)

    def test_validate_threads_bad(self):
        with self.assertRaises(SystemExit):
            validate_threads("o")


class TestParseMmseqsVersion(unittest.TestCase):
    """Tests of parse_mmseqs_version (issue #451)"""

    def test_parse_known_release(self):
        self.assertEqual(parse_mmseqs_version("18.8cc5c"), (18, "8cc5c"))

    def test_parse_newer_release(self):
        self.assertEqual(parse_mmseqs_version("19.ab57c"), (19, "ab57c"))

    def test_parse_old_release(self):
        self.assertEqual(parse_mmseqs_version("13.45111"), (13, "45111"))

    def test_parse_hash_only_build(self):
        self.assertEqual(
            parse_mmseqs_version("7e2840992948ee89dcc336522dc98a74fe0adf00"),
            (14, "7e284"),
        )

    def test_parse_unknown(self):
        self.assertIsNone(parse_mmseqs_version("abcdef0123456789"))


def _fake_mmseqs_popen(version):
    """sp.Popen that reports the given MMseqs2 version, otherwise unchanged"""
    real_popen = sp.Popen

    def fake_popen(cmd, *args, **kwargs):
        if cmd[:2] == ["mmseqs", "version"]:
            cmd = ["echo", version]
        return real_popen(cmd, *args, **kwargs)

    return fake_popen


class TestCheckDependenciesMmseqs(unittest.TestCase):
    """check_dependencies must not crash on new or unrecognised MMseqs2 versions (issue #451)"""

    def test_mmseqs_newer_release(self):
        with mock.patch(
            "pharokka.input_commands.sp.Popen", _fake_mmseqs_popen("19.ab57c")
        ):
            check_dependencies(skip_mash=True)

    def test_mmseqs_unrecognised_version(self):
        with mock.patch(
            "pharokka.input_commands.sp.Popen",
            _fake_mmseqs_popen("abcdef0123456789"),
        ):
            check_dependencies(skip_mash=True)
