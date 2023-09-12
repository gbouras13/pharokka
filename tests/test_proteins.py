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

from bin.util import remove_directory

# import functions


# test data
test_data = Path("tests/test_data")
database_dir = Path(f"{test_data}/database")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{test_data}/Meta_example")
proteins_data = Path(f"{test_data}/proteins")
threads = 4

logger.add(lambda _: sys.exit(1), level="ERROR")


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


temp_dir = Path(f"{proteins_data}/fake_out")


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


def test_download(tmp_dir):
    """test pharokka download"""
    cmd = f"install_databases.py -o {database_dir}"
    exec_command(cmd)


def test_proteins(tmp_dir):
    """test pharokka proteins"""
    input_fasta: Path = f"{proteins_data}/phanotate.faa"
    cmd = f"pharokka_proteins.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f"
    exec_command(cmd)


def test_proteins_hmm_only(tmp_dir):
    """test pharokka proteins hmm_only"""
    input_fasta: Path = f"{proteins_data}/phanotate.faa"
    cmd = f"pharokka_proteins.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --hmm_only"
    exec_command(cmd)


def test_proteins_mmseqs_only(tmp_dir):
    """test pharokka proteins mmseqs_only"""
    input_fasta: Path = f"{proteins_data}/phanotate.faa"
    cmd = f"pharokka_proteins.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --mmseqs2_only"
    exec_command(cmd)


temp_dir = Path(f"{test_data}/fake_out")


class testFails(unittest.TestCase):
    """Tests for fails"""

    def test_proteins_pred(self):
        """tests that pharokka breaks if nucleotide input"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{meta_data}/fake_meta.fasta"
            cmd = f"pharokka_proteins.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t {threads} -f"
            exec_command(cmd)


remove_directory(f"{temp_dir}")
remove_directory(f"{database_dir}")
