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

# test_data = Path("tests/test_data")
# overall_test_data = Path(f"{test_data}/overall_inputs")


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


# def test_phage(tmp_dir):
#     """test phage"""
#     input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#     cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_chrom(tmp_dir):
#     """test chrom"""
#     input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
#     cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_chrom_diff_eval(tmp_dir):
#     """test chrom with different e value"""
#     input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
#     cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -e 0.1 -f"
#     exec_command(cmd)


# def test_chrom_mystery_autocomplete(tmp_dir):
#     """test chrom - with phage (so no hit)"""
#     input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#     cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
#     exec_command(cmd)


# def test_chrom_nearest_autocomplete(tmp_dir):
#     """test chrom - with phage (so no hit)"""
#     input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#     cmd = f"dnaapler chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
#     exec_command(cmd)


# def test_plas(tmp_dir):
#     """test plas"""
#     input_fasta: Path = f"{overall_test_data}/plasmid.fasta"
#     cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_plas_nearest_autocomplete(tmp_dir):
#     """test plas - with phage (so no hit)"""
#     input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#     cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
#     exec_command(cmd)


# def test_plas_mystery_autocomplete(tmp_dir):
#     """test plas - with phage (so no hit)"""
#     input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#     cmd = f"dnaapler plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
#     exec_command(cmd)


# def test_phage_mystery(tmp_dir):
#     """test phage mystery - with plasmid (so no hit)"""
#     input_fasta: Path = f"{overall_test_data}/plasmid.fasta"
#     cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f -a mystery"
#     exec_command(cmd)


# def test_phage_nearest(tmp_dir):
#     """test phage nearest - with plasmid (so no hit)"""
#     input_fasta: Path = f"{overall_test_data}/plasmid.fasta"
#     cmd = f"dnaapler phage -i {input_fasta} -o {tmp_dir} -t 1 -f -a nearest"
#     exec_command(cmd)


# def test_mys(tmp_dir):
#     """test mystery"""
#     input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
#     cmd = f"dnaapler mystery -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_nearest(tmp_dir):
#     """test nearest"""
#     input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
#     cmd = f"dnaapler nearest -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_custom(tmp_dir):
#     """test custom"""
#     input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
#     custom: Path = f"{test_data}/fake_custom.faa"
#     cmd = f"dnaapler custom -i {input_fasta} -o {tmp_dir} -c {custom} -t 1 -f"
#     exec_command(cmd)


# def test_bulk_chromosome(tmp_dir):
#     """test bulk chromosome"""
#     input_fasta: Path = f"{overall_test_data}/bulk_chromosome.fasta"
#     cmd = f"dnaapler bulk -m chromosome -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_bulk_phage(tmp_dir):
#     """test bulk phage"""
#     input_fasta: Path = f"{overall_test_data}/bulk_phage.fasta"
#     cmd = f"dnaapler bulk -m phage -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_bulk_plasmid(tmp_dir):
#     """test bulk plasmid"""
#     input_fasta: Path = f"{overall_test_data}/bulk_plasmid.fasta"
#     cmd = f"dnaapler bulk -m plasmid -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_bulk_custom(tmp_dir):
#     """test bulk custom"""
#     input_fasta: Path = f"{overall_test_data}/bulk_chromosome.fasta"
#     custom: Path = f"{test_data}/fake_custom.faa"
#     cmd = f"dnaapler bulk -m custom -i {input_fasta} -o {tmp_dir} -c {custom} -t 1 -f"
#     exec_command(cmd)


# def test_all(tmp_dir):
#     """test all"""
#     input_fasta: Path = f"{overall_test_data}/all_test.fasta"
#     cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f"
#     exec_command(cmd)


# def test_all_ignore(tmp_dir):
#     """test all"""
#     input_fasta: Path = f"{overall_test_data}/all_test.fasta"
#     ignore_file: Path = f"{overall_test_data}/ignore.txt"
#     cmd = f"dnaapler all  -i {input_fasta} -o {tmp_dir} -t 1 -f --ignore {ignore_file}"
#     exec_command(cmd)


# def test_citation():
#     """test citation"""
#     cmd = f"dnaapler citation"
#     exec_command(cmd)


# class TestExits(unittest.TestCase):
#     """Tests of End to End common failures"""

#     def test_chrom_double(self):
#         """test chrom double to test force is working"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/chromosome.fasta"
#             outdir: Path = f"{overall_test_data}/chrom_out"
#             cmd = f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1"
#             exec_command(cmd)
#             exec_command(cmd)

#     def test_phage_chrom(self):
#         """test chrom with phage input - test where nothing found"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#             outdir: Path = f"{overall_test_data}/phage_out"
#             cmd = f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1 -f"
#             exec_command(cmd)

#     def test_chrom_bad_autocomplete(self):
#         """test chrom with bad autcomplete"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/SAOMS1.fasta"
#             outdir: Path = f"{overall_test_data}/phage_out"
#             cmd = (
#                 f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1 -f -a bad_auto"
#             )
#             exec_command(cmd)

#     def test_chrom_already_oriented(self):
#         """test chrom with already oriented input"""
#         with self.assertRaises(RuntimeError):
#             input_fasta: Path = f"{overall_test_data}/SAOMS1_reoriented.fasta"
#             outdir: Path = f"{overall_test_data}/chrom_out"
#             cmd = f"dnaapler chromosome -i {input_fasta} -o {outdir} -t 1 -f"
#             exec_command(cmd)

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


# remove_directory(f"{overall_test_data}/phage_out")
# remove_directory(f"{overall_test_data}/chrom_out")
# remove_directory(f"{overall_test_data}/bulk_out")
