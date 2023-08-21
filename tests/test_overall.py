# """
# Unit tests for pharokka overall

# Usage: pytest .

# """

# # import
# import os
# import shutil

# # import functions
# import subprocess
# import unittest
# from pathlib import Path

# import pytest


# import sys
# # import
# import unittest
# from pathlib import Path
# from unittest.mock import patch
# from loguru import logger

# from lib.input_commands import (instantiate_dirs, validate_fasta,
#                                 validate_gene_predictor, validate_meta,
#                                 validate_strand, validate_terminase,
#                                 validate_terminase_start, validate_threads)
# from lib.util import remove_directory

# # import functions


# # test data
# test_data = Path("tests/test_data")
# functions_data = Path(f"{test_data}/functions_files")
# database_dir = Path(f"{test_data}/functions_data")
# overall_data = Path(f"{test_data}/overall")
# meta_data = Path(f"{overall_data}/Meta_example")
# standard_data = Path(f"{overall_data}/Standard_examples")
# logger.add(lambda _: sys.exit(1), level="ERROR")


# def remove_directory(dir_path):
#     if os.path.exists(dir_path):
#         shutil.rmtree(dir_path)


# @pytest.fixture(scope="session")
# def tmp_dir(tmpdir_factory):
#     return tmpdir_factory.mktemp("tmp")


# def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
#     """executes shell command and returns stdout if completes exit code 0
#     Parameters
#     ----------
#     cmnd : str
#       shell command to be executed
#     stdout, stderr : streams
#       Default value (PIPE) intercepts process output, setting to None
#       blocks this."""

#     proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
#     out, err = proc.communicate()
#     if proc.returncode != 0:
#         raise RuntimeError(f"FAILED: {cmnd}\n{err}")
#     return out.decode("utf8") if out is not None else None

# def test_download(tmp_dir):
#     """test pharokka download"""
#     cmd = f"install_databases.py -o {database_dir}"
#     exec_command(cmd)

# def test_overall(tmp_dir):
#     """test pharokka overall"""
#     input_fasta: Path = f"{standard_data}/SAOMS1.fasta"
#     cmd = f"pharokka.py phage -i {input_fasta} -d {database_dir} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


#     def test_bulk_single_genome(self):
#         """test bulk with single genome"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/SAOMS1_reoriented.fasta"
#             outdir: Path = f"{overall_test_data}/chrom_out"
#             cmd = f"dnaapler bulk -m phage -i {input_fasta} -o {outdir} -t 1 -f"
#             exec_command(cmd)

#     def test_bulk_bad_mode(self):
#         """test bulk with incorrect -m value"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/SAOMS1_reoriented.fasta"
#             outdir: Path = f"{overall_test_data}/chrom_out"
#             cmd = f"dnaapler bulk -m bad -i {input_fasta} -o {outdir} -t 1 -f"
#             exec_command(cmd)

#     def test_bulk_custom_no_db(self):
#         """test bulk custom with no db"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/bulk_chromosome.fasta"
#             outdir: Path = f"{overall_test_data}/chrom_out"
#             cmd = f"dnaapler bulk -m custom -i {input_fasta} -o {outdir} -t 1 -f"
#             exec_command(cmd)

#     def test_bulk_phage_with_chrom(self):
#         """test bulk chromosome with phage genomes"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/bulk_phage.fasta"
#             outdir: Path = f"{overall_test_data}/chrom_out"
#             cmd = f"dnaapler bulk -m chromosome -i {input_fasta} -o {outdir} -t 1 -f"
#             exec_command(cmd)

#     def test_all_no_hits(self):
#         """test all with no blast hits"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{test_data}/nucl_test.fna"
#             outdir: Path = f"{overall_test_data}/bulk_out"
#             cmd = f"dnaapler all -i {input_fasta} -o {outdir} -t 1 -f "
#             exec_command(cmd)


# remove_directory(f"{database_dir}")
# remove_directory(f"{overall_test_data}/chrom_out")
# remove_directory(f"{overall_test_data}/bulk_out")
