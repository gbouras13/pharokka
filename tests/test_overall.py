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
functions_data = Path(f"{test_data}/functions_files")
database_dir = Path(f"{test_data}/database")
overall_data = Path(f"{test_data}/overall")
meta_data = Path(f"{overall_data}/Meta_example")
standard_data = Path(f"{overall_data}/Standard_examples")
standard_data_output = Path(f"{standard_data}/SAOMS1_Output")
stop_recoding_data = Path(f"{overall_data}/stop_recoding")
custom_db = Path(f"{test_data}/custom_db/microvirus.h3m")
custom_data = Path(f"{overall_data}/custom_examples")
tmrna_data = Path(f"{overall_data}/tmRNA_example")
AMR_data = Path(f"{overall_data}/AMR_example")
CRISPR_data = Path(f"{overall_data}/CRISPR_example")
VFDB_data = Path(f"{overall_data}/VFDB_example")
genbank_data = Path(f"{overall_data}/genbank_examples")
logger.add(lambda _: sys.exit(1), level="ERROR")
threads = 8


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


temp_dir = Path(f"{test_data}/fake_out")


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


def test_overall(tmp_dir):
    """test pharokka overall"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f"
    exec_command(cmd)

def test_overall_mash_distance(tmp_dir):
    """test pharokka overall with stricter mash distance"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --mash_distance 0.05"
    exec_command(cmd)


def test_overall_crispr(tmp_dir):
    """test pharokka overall crispr"""
    input_fasta: Path = f"{CRISPR_data}/Biggiephage_A_fullcontig_CasΦ1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f"
    exec_command(cmd)


def test_overall_crispr_minced_args(tmp_dir):
    """test pharokka crispr with minced args"""
    input_fasta: Path = f"{CRISPR_data}/Biggiephage_A_fullcontig_CasΦ1.fasta"
    cmd = f'pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -g prodigal --minced_args "minNR 2 -minRL 21" '
    exec_command(cmd)


def test_overall_vfdb(tmp_dir):
    """test pharokka overall on a phage with vfdb hits. Also include --skip_extra_annotations"""
    input_fasta: Path = f"{VFDB_data}/NC_004617.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --skip_extra_annotations"
    exec_command(cmd)


def test_overall_amr(tmp_dir):
    """test pharokka overall amr also includes '--skip_mash"""
    input_fasta: Path = f"{AMR_data}/NC_007458.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --skip_mash"
    exec_command(cmd)


def test_overall_tmrna(tmp_dir):
    """test pharokka overall tmrna"""
    input_fasta: Path = f"{tmrna_data}/NC_051700.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f"
    exec_command(cmd)


def test_overall_numeric_header_prodigal(tmp_dir):
    """#317 test pharokka prodigal with  numeric header #334"""
    input_fasta: Path = f"{standard_data}/SAOMS1_numeric_header.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -g prodigal --fast -m"
    exec_command(cmd)


def test_overall_numeric_header_phanotate(tmp_dir):
    """#317 test pharokka phanotate  with numeric header #334"""
    input_fasta: Path = f"{standard_data}/SAOMS1_numeric_header.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --fast "
    exec_command(cmd)


def test_overall_numeric_header_prodigal_gv(tmp_dir):
    """#317 test pharokka prodigal-gv with numeric header #334"""
    input_fasta: Path = f"{standard_data}/SAOMS1_numeric_header.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m --fast"
    exec_command(cmd)


def test_meta(tmp_dir):
    """test pharokka meta"""
    input_fasta: Path = f"{meta_data}/combined_meta.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m"
    exec_command(cmd)


def test_meta_unicycler_header_prodigal(tmp_dir):
    """#317 test pharokka overall with unicycler header - with all trna, tmras crispr etc"""
    input_fasta: Path = f"{meta_data}/combined_meta_unicycler_headers.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -g prodigal-gv -m"
    exec_command(cmd)


def test_meta_prodigal_gv(tmp_dir):
    """test pharokka meta with prodigal-gv"""
    input_fasta: Path = f"{meta_data}/combined_meta.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m -g prodigal-gv"
    exec_command(cmd)


def test_meta_dnaapler_all_bug(tmp_dir):
    """test pharokka meta dnaapler bug and split"""
    input_fasta: Path = f"{meta_data}/combined_meta.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m -s --dnaapler --meta_hmm"
    exec_command(cmd)


def test_overall_locus(tmp_dir):
    """test pharokka overall locus tag prefix"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -l SAOMS1 -p SAOMS1"
    exec_command(cmd)


def test_custom(tmp_dir):
    """test pharokka overall with custom db"""
    input_fasta: Path = f"{custom_data}/MH649026.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} --custom_hmm {custom_db} -o {tmp_dir} -t {threads} -f "
    exec_command(cmd)


def test_custom_meta(tmp_dir):
    """test pharokka overall with custom db"""
    input_fasta: Path = f"{custom_data}/hundred_microviruses.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} --custom_hmm {custom_db} -o {tmp_dir} -t {threads} -f -m"
    exec_command(cmd)


def test_overall_prodigal(tmp_dir):
    """test pharokka overall prodigal"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -g prodigal"
    exec_command(cmd)


def test_overall_stop_recode(tmp_dir):
    """test pharokka overall recoded"""
    input_fasta: Path = f"{stop_recoding_data}/table_4/SRR1747055_scaffold_7.fa"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -g prodigal -c 4"
    exec_command(cmd)


def test_overall_dnaapler(tmp_dir):
    """test pharokka overall dnaapler"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --dnaapler"
    exec_command(cmd)


def test_overall_fast(tmp_dir):
    """test pharokka overall fast"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --fast"
    exec_command(cmd)


def test_overall_mmseqs(tmp_dir):
    """test pharokka overall mmseqs2_only"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --mmseqs2_only"
    exec_command(cmd)


def test_meta_no_cds_contig(tmp_dir):
    """test pharokka meta with a contig with empty contigs"""
    input_fasta: Path = f"{meta_data}/fake_meta.fa"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m"
    exec_command(cmd)


######
# pharokka CI was timing out (>6 hours)
# These are covered by other rules anyway
# so just run as is

# def test_meta_hmm(tmp_dir):
#     """test pharokka meta hmm"""
#     input_fasta: Path = f"{meta_data}/fake_meta.fa"
#     cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m --meta_hmm"
#     exec_command(cmd)


# def test_meta_split(tmp_dir):
#     """test pharokka meta split"""
#     input_fasta: Path = f"{meta_data}/fake_meta.fa"
#     cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f -m -s"
#     exec_command(cmd)


def test_terminase(tmp_dir):
    """test pharokka terminase"""
    input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
    cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {tmp_dir} -t {threads} -f --terminase --terminase_start 340 --terminase_strand neg"
    exec_command(cmd)


def test_overall_genbank(tmp_dir):
    """test pharokka overall with genbank input"""
    input_gbk: Path = f"{genbank_data}/SAOMS1.gbk"
    cmd = f"pharokka.py -i {input_gbk} -d {database_dir} -o {tmp_dir} -t {threads} -f --genbank"
    exec_command(cmd)


def test_overall_genbank_meta(tmp_dir):
    """test pharokka overall meta with genbank input"""
    input_gbk: Path = f"{genbank_data}/hundred_microviruses.gbk"
    cmd = f"pharokka.py -i {input_gbk} -d {database_dir} -o {tmp_dir} -t {threads} -f --genbank -m --meta_hmm"
    exec_command(cmd)


class testFails(unittest.TestCase):
    """Tests for fails"""

    def test_dupe_header(self):
        """tests that pharokka exits if a duplicate header is passed"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/dupe_header.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m"
            exec_command(cmd)

    def test_overall_hash_header(self):
        """test pharokka overall with # in header - should error out"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1_hash_header.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t {threads} -f"
            exec_command(cmd)

    def test_meta_with_single_contig(self):
        """tests that pharokka exits if single contig is passed to meta"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m"
            exec_command(cmd)

    def test_meta_terminas(self):
        """tests that pharokka exits if multiple contigs passed to meta"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{meta_data}/fake_meta.fa"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f -m --terminase"
            exec_command(cmd)

    def test_bad_terminase(self):
        """tests that pharokka exits if only terminase specified"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f --terminase"
            exec_command(cmd)

    def test_bad_terminase_start(self):
        """tests that pharokka exits if bad terminase start specified"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f --terminase --terminase_start sf"
            exec_command(cmd)

    def test_bad_terminase_strand(self):
        """tests that pharokka exits if bad strand"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -f --terminase --terminase_start 34 --terminase_strand posit"
            exec_command(cmd)

    def test_bad_threads(self):
        """tests that pharokka exits if bad threads"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = (
                f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t sf -f"
            )
            exec_command(cmd)

    def test_bad_fast_mmseqs2_only(self):
        """tests that pharokka exits if both fast and mmseqs2_only"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t {threads} --fast --mmseqs2_only -f"
            exec_command(cmd)

    def test_bad_fast_meta_gmm(self):
        """tests that pharokka exits if both fast and meta_hmm"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t {threads} --fast --meta_hmm -m -f"
            exec_command(cmd)

    def test_bad_fast_meta_gmm(self):
        """tests that pharokka exits if both mmseqs2_only and meta_hmm"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t {threads} --mmseqs2_only --meta_hmm -m -f"
            exec_command(cmd)

    def test_bad_gene_pred(self):
        """tests that pharokka exits if bad gene predictior"""
        with self.assertRaises(RuntimeError):
            input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} -o {temp_dir} -t 1 -g proddigal -f"
            exec_command(cmd)

    def test_bad_genbank(self):
        """test pharokka overall with genbank input but no --genbank"""
        with self.assertRaises(RuntimeError):
            input_gbk: Path = f"{standard_data_output}/SAOMS1.fasta"
            cmd = f"pharokka.py -i {input_gbk} -d {database_dir} -o {temp_dir} -t {threads} -f "
            exec_command(cmd)

    def test_custom_bad(self):
        """test pharokka overall with bad custom db"""
        with self.assertRaises(RuntimeError):
            custom_db: Path = f"{standard_data_output}/SAOMS1.fasta"
            input_fasta: Path = f"{custom_data}/MH649026.fasta"
            cmd = f"pharokka.py -i {input_fasta} -d {database_dir} --custom_hmm {custom_db} -o {temp_dir} -t {threads} -f "
            exec_command(cmd)


remove_directory(temp_dir)
remove_directory(f"{database_dir}")
