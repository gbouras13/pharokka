"""
Unit tests for plassembler.

Usage: pytest

"""


# import
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest
from loguru import logger

from bin.processes import (run_aragorn, run_mash_sketch, run_minced,
                           run_phanotate, run_pyrodigal)
# import functions
from bin.util import remove_directory

# test data
test_data = Path("tests/test_data")
functions_data = Path(f"{test_data}/functions_files")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{overall_data}/Meta_example")
standard_data = Path(f"{overall_data}/Standard_examples")
standard_data_output = Path(f"{standard_data}/SAOMS1_Output")
logdir = Path(f"{test_data}/logs")


logger.add(lambda _: sys.exit(1), level="ERROR")


# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class testGenePred(unittest.TestCase):
    """Tests for gene predictors"""

    def test_run_phanotate(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        run_phanotate(fasta, standard_data_output, logdir)

    def test_run_pyrodigal(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 11
        meta = False
        run_pyrodigal(fasta, standard_data_output, meta, coding_table)  # meta = False

    def test_run_pyrodigal_meta(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 11
        meta = True
        run_pyrodigal(fasta, standard_data_output, meta, coding_table)  # meta = False

    def test_run_pyrodigal_c4(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 4
        meta = False
        run_pyrodigal(fasta, standard_data_output, meta, coding_table)  # meta = False

    def test_run_minced(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        prefix = "pharokka"
        run_minced(fasta, standard_data_output, prefix, logdir)

    def test_run_aragorn(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        prefix = "pharokka"
        run_aragorn(fasta, standard_data_output, prefix, logdir)

    def test_run_aragorn(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        run_mash_sketch(fasta, standard_data_output, logdir)


remove_directory(logdir)
