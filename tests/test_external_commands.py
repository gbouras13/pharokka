"""
Unit tests for plassembler.

Usage: pytest

"""

# import
import sys
from pathlib import Path

import pytest
from loguru import logger

from pharokka.processes import (
    run_aragorn,
    run_mash_sketch,
    run_minced,
    run_phanotate,
    run_pyrodigal,
    run_pyrodigal_gv,
)

# import functions
from pharokka.util import remove_directory

# test data
test_data = Path("tests/test_data")
functions_data = Path(f"{test_data}/functions_files")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{overall_data}/Meta_example")
standard_data = Path(f"{overall_data}/Standard_examples")
logdir = Path(f"{test_data}/logs")


logger.add(lambda _: sys.exit(1), level="ERROR")


# make fake tempdir for testing
@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class TestGenePred:
    """Tests for gene predictors"""

    def test_run_phanotate(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        run_phanotate(fasta, tmp_dir, logdir)

    def test_run_pyrodigal(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 11
        meta = False
        threads = 2
        run_pyrodigal(fasta, tmp_dir, meta, coding_table, threads)  # meta = False

    def test_run_pyrodigal_small_fasta(self, tmp_dir):
        """
        handle instance where input is under 20000 nucleotides
        """
        fasta: Path = f"{standard_data}/SAOMS1_subset.fasta"
        coding_table = 11
        meta = False
        threads = 2
        run_pyrodigal(fasta, tmp_dir, meta, coding_table, threads)  # meta = False

    def test_run_pyrodigal_meta(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 11
        meta = True
        threads = 2
        run_pyrodigal(fasta, tmp_dir, meta, coding_table, threads)  # meta = False

    def test_run_pyrodigal_c4(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        coding_table = 4
        meta = False
        threads = 2
        run_pyrodigal(fasta, tmp_dir, meta, coding_table, threads)  # meta = False

    def test_run_pyrodigal_gv(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        threads = 2
        run_pyrodigal_gv(fasta, tmp_dir, threads)  # meta = False

    def test_run_minced(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        prefix = "pharokka"
        minced_args = "minNR 2 -minRL 21"
        run_minced(fasta, tmp_dir, prefix, minced_args, logdir)

    def test_run_aragorn(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        prefix = "pharokka"
        run_aragorn(fasta, tmp_dir, prefix, logdir)

    def test_run_mash_sketch(self, tmp_dir):
        fasta: Path = f"{standard_data}/SAOMS1.fasta"
        run_mash_sketch(fasta, tmp_dir, logdir)


remove_directory(logdir)
