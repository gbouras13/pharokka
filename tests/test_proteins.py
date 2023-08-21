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
from pathlib import Path
from unittest.mock import patch
import unittest
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

database_dir = Path(f"{test_data}/database")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{overall_data}/Meta_example")
standard_data = Path(f"{overall_data}/Standard_examples")

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


def test_proteins(tmp_dir):
    """test pharokka overall"""
    input_fasta: Path = f"{meta_data}/fake_meta.fa"
    cmd = (
        f"pharokka_proteins.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f"
    )
    exec_command(cmd)

temp_dir = Path(f"{test_data}/fake_out")

class testFails(unittest.TestCase):
    """Tests for fails"""
    def test_proteins_pred(self):
        """tests that pharokka exits if bad gene predictior"""
        with self.assertRaises(SystemExit):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka_proteins.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f"
            exec_command(cmd)


remove_directory(f"{temp_dir}")
# remove_directory(f"{database_dir}")