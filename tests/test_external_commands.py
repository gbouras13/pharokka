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
                           run_phanotate, run_pyrodigal, run_pyrodigal_gv)
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
        threads = 2
        run_pyrodigal(
            fasta, standard_data_output, meta, coding_table, threads
        )  # meta = False

    def test_run_pyrodigal_small_fasta(self):
        """
        handle instance where input is under 20000 nucleotides
        """
        fasta: Path = f"{standard_data}/SAOMS1_subset.fasta"
        coding_table = 11
        meta = False
        threads = 2
        run_pyrodigal(
            fasta, standard_data_output, meta, coding_table, threads
        )  # meta = False

    def test_run_pyrodigal_meta(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 11
        meta = True
        threads = 2
        run_pyrodigal(
            fasta, standard_data_output, meta, coding_table, threads
        )  # meta = False

    def test_run_pyrodigal_c4(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 4
        meta = False
        threads = 2
        run_pyrodigal(
            fasta, standard_data_output, meta, coding_table, threads
        )  # meta = False

    def test_run_pyrodigal_gv(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        threads = 2
        run_pyrodigal_gv(fasta, standard_data_output, threads)  # meta = False

    def test_run_minced(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        prefix = "pharokka"
        minced_args = "minNR 2 -minRL 21"
        run_minced(fasta, standard_data_output, prefix, minced_args, logdir)

    def test_run_aragorn(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        prefix = "pharokka"
        run_aragorn(fasta, standard_data_output, prefix, logdir)

    def test_run_aragorn(self):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        run_mash_sketch(fasta, standard_data_output, logdir)


remove_directory(logdir)
