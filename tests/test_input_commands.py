"""
Unit tests for pharokka

Usage: pytest

"""

import sys
# import
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest
from loguru import logger

from bin.input_commands import (instantiate_dirs, validate_fasta,
                                validate_gene_predictor, validate_meta,
                                validate_strand, validate_terminase,
                                validate_terminase_start, validate_threads)
from bin.util import remove_directory

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


class TestValidateTerminase(unittest.TestCase):
    """Tests of validate_threads"""

    def test_validate_threads(self):
        validate_threads(1)

    def test_validate_threads_bad(self):
        with self.assertRaises(SystemExit):
            validate_threads("o")
