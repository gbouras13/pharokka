"""
Unit tests for pharokka overall

Usage: pytest .

"""

# import
import os
import shutil
# import functions
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import pytest
from loguru import logger

from lib.input_commands import (instantiate_dirs, validate_fasta,
                                validate_gene_predictor, validate_meta,
                                validate_strand, validate_terminase,
                                validate_terminase_start, validate_threads)
from lib.util import remove_directory

# import functions


# test data
test_data = Path("tests/test_data")
functions_data = Path(f"{test_data}/functions_files")
database_dir = Path(f"{test_data}/database")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{overall_data}/Meta_example")
standard_data = Path(f"{overall_data}/Standard_examples")
stop_recoding_data = Path(f"{overall_data}/stop_recoding")
logger.add(lambda _: sys.exit(1), level="ERROR")


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


# def test_download(tmp_dir):
#     """test pharokka download"""
#     cmd = f"install_databases.py -o {database_dir}"
#     exec_command(cmd)


def test_overall(tmp_dir):
    """test pharokka overall"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f"
    exec_command(cmd)


def test_overall_prodigal(tmp_dir):
    """test pharokka overall prodigal"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f -g prodigal"
    exec_command(cmd)


def test_overall_stop_recode(tmp_dir):
    """test pharokka overall recoded"""
    input_fasta: Path = f"{stop_recoding_data}/table_4/SRR1747055_scaffold_7.fa"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f -g prodigal -c 4"
    exec_command(cmd)


def test_overall_dnaapler(tmp_dir):
    """test pharokka overall dnaapler"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f --dnaapler"
    exec_command(cmd)


def test_overall_fast(tmp_dir):
    """test pharokka overall fast"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f --fast"
    exec_command(cmd)


def test_overall_mmseqs(tmp_dir):
    """test pharokka overall mmseqs2_only"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f --mmseqs2_only"
    exec_command(cmd)


def test_meta(tmp_dir):
    """test pharokka meta"""
    input_fasta: Path = f"{meta_data}/fake_meta.fa"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f -m"
    exec_command(cmd)


def test_meta_split(tmp_dir):
    """test pharokka meta"""
    input_fasta: Path = f"{meta_data}/fake_meta.fa"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f -m -s"
    exec_command(cmd)


def test_terminase(tmp_dir):
    """test pharokka terminase"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f --terminase --terminase_start 340 --terminase_strand neg"
    exec_command(cmd)



temp_dir = Path(f"{test_data}/fake_out")

class testFails(unittest.TestCase):
    """Tests for fails"""
    

    def test_meta_with_single_contig(self):
        """tests that pharokka exits if single contig is passed to meta"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m"
            exec_command(cmd)

    def test_meta_terminas(self):
        """tests that pharokka exits if multiple contigs passed to meta"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{meta_data}/fake_meta.fa"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m --terminase"
            exec_command(cmd)

    def test_bad_terminase(self):
        """tests that pharokka exits if only terminase specified"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f --terminase"
            exec_command(cmd)

    def test_bad_terminase_start(self):
        """tests that pharokka exits if bad terminase start specified"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f --terminase --terminase_start sf"
            exec_command(cmd)

    def test_bad_terminase_strand(self):
        """tests that pharokka exits if bad strand"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f --terminase --terminase_start 34 --terminase_strand posit"
            exec_command(cmd)

    def test_bad_threads(self):
        """tests that pharokka exits if bad threads"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = (
                f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t sf -f"
            )
            exec_command(cmd)

    def test_bad_gene_pred(self):
        """tests that pharokka exits if bad gene predictior"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -g proddigal -f"
            exec_command(cmd)


remove_directory(temp_dir)
# remove_directory(f"{database_dir}")
