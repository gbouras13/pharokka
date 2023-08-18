"""
Unit tests for pharokka

Usage: pytest

"""

import glob
import os
import sys

# import
import unittest
from pathlib import Path
from unittest.mock import patch

import click
import pandas as pd
import pytest
from loguru import logger

from lib.input_commands import ( instantiate_dirs, validate_fasta)
from lib.util import remove_directory

# import functions


# test data
test_data = Path("tests/test_data")
functions_data = Path(f"{test_data}/functions_files")
logger.add(lambda _: sys.exit(1), level="ERROR")

@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


class TestInstantiateDir(unittest.TestCase):
    """Tests of instantiate dir"""

    def test_instantiate_dirs(self):
        fake_out_dir: Path = f"{test_data}/fake_out"
        instantiate_dirs(fake_out_dir, meta = False, force = False)
        remove_directory(fake_out_dir)
        
    def test_instantiate_overwrite(self):
        fake_out_dir_for_test: Path = f"{test_data}/fake_out_for_test"
        instantiate_dirs(fake_out_dir_for_test, meta = False, force = True)

    def test_instantiate_overwrite_fail(self):
        with self.assertRaises(SystemExit):
            fake_out_dir_for_test: Path = f"{test_data}/fake_out_for_test"
            instantiate_dirs(fake_out_dir_for_test, meta = False, force = False)

    def test_instantiate_dirs_fail_output_exists(self):
        with self.assertRaises(SystemExit):
            instantiate_dirs(test_data, meta = False, force = False)

class TestValidateFasta(unittest.TestCase):
    """Tests of instantiate dir"""

    def test_validate_fasta(self):
        fasta: Path = f"{functions_data}/good.fasta"
        validate_fasta(fasta)
        
    def test_validate_fasta_bad(self):
        with self.assertRaises(SystemExit):
            fasta: Path = f"{functions_data}/bad.fasta"
            validate_fasta(fasta)

    


        # checks the ctx is the same, no error

#     def test_non_fasta_custom_input(self):
#         with self.assertRaises(SystemExit):
#             nucleotide_fasta_file = os.path.join(test_data, "non_fasta.txt")
#             validate_custom_db_fasta(nucleotide_fasta_file)

#     def test_non_aa_custom_fasta_input(self):
#         with self.assertRaises(SystemExit):
#             nucleotide_fasta_file = os.path.join(test_data, "nucl_test.fna")
#             validate_custom_db_fasta(nucleotide_fasta_file)


# class TestReorientSequence(unittest.TestCase):
#     """Tests for reorient_sequence"""

#     def test_reorient_sequence_outside_range(self):
#         # Test scenario where the row is outside of range
#         with self.assertRaises(KeyError):
#             blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
#             col_list = [
#                 "qseqid",
#                 "qlen",
#                 "sseqid",
#                 "slen",
#                 "length",
#                 "qstart",
#                 "qend",
#                 "sstart",
#                 "send",
#                 "pident",
#                 "nident",
#                 "gaps",
#                 "mismatch",
#                 "evalue",
#                 "bitscore",
#                 "qseq",
#                 "sseq",
#             ]
#             blast_df = pd.read_csv(
#                 blast_file, delimiter="\t", index_col=False, names=col_list
#             )
#             input = os.path.join(test_data, "SAOMS1.fasta")
#             out_file = os.path.join(test_data, "fake_reoriented.fasta")
#             gene = "terL"
#             i = int(3)
#             reorient_sequence(blast_df, input, out_file, gene, i)

#     def test_reorient_sequence_correct(self):
#         # test where it works as expected
#         blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
#         col_list = [
#             "qseqid",
#             "qlen",
#             "sseqid",
#             "slen",
#             "length",
#             "qstart",
#             "qend",
#             "sstart",
#             "send",
#             "pident",
#             "nident",
#             "gaps",
#             "mismatch",
#             "evalue",
#             "bitscore",
#             "qseq",
#             "sseq",
#         ]
#         blast_df = pd.read_csv(
#             blast_file, delimiter="\t", index_col=False, names=col_list
#         )
#         input = os.path.join(test_data, "SAOMS1.fasta")
#         out_file = os.path.join(test_data, "fake_reoriented.fasta")
#         gene = "terL"
#         i = int(1)
#         reorient_sequence(blast_df, input, out_file, gene, i)


# class TestReorientSequenceRandom(unittest.TestCase):
#     """Tests for reorient_sequence_random"""

#     def test_reorient_sequence_random_bad_strand(self):
#         # Test scenario where the strand is outside 1 or -1
#         with self.assertRaises(UnboundLocalError):
#             input = os.path.join(test_data, "SAOMS1.fasta")
#             out_file = os.path.join(test_data, "fake_reoriented.fasta")
#             start = 1000
#             strand = 24
#             reorient_sequence_random(input, out_file, start, strand)


# class TestBlastOutput(unittest.TestCase):
#     """Tests for process_blast_output_and_reorient"""

#     def test_process_blast_output_and_reorient_invalid_blast_file(self):
#         # Test scenario where the blast input is invalud
#         with self.assertRaises(SystemExit):
#             blast_file = pd.DataFrame({"qstart": [1]})
#             input = os.path.join(test_data, "SAOMS1.fasta")
#             output = os.path.join(test_data, "fake_reoriented.fasta")
#             gene = "terL"
#             process_blast_output_and_reorient(input, blast_file, output, gene)

#     def test_process_blast_output_and_reorient_already_oriented(self):
#         # Test scenario where the blast output suggests the contig is already oriented correctly
#         blast_file = os.path.join(test_data, "SAOMS1_blast_output_already_oriented.txt")
#         input = os.path.join(test_data, "SAOMS1.fasta")
#         output = os.path.join(test_data, "fake_reoriented.fasta")
#         gene = "terL"
#         blast_success = process_blast_output_and_reorient(
#             input, blast_file, output, gene
#         )
#         assert blast_success == True

#     def test_process_blast_output_and_reorient_wrong_start_codon(self):
#         # Test scenario where the best BLAST hit has no valid start codon
#         blast_file = os.path.join(
#             test_data, "SAOMS1_blast_output_wrong_start_codon.txt"
#         )
#         input = os.path.join(test_data, "SAOMS1.fasta")
#         output = os.path.join(test_data, "fake_reoriented.fasta")
#         gene = "terL"
#         blast_success = process_blast_output_and_reorient(
#             input, blast_file, output, gene
#         )
#         assert blast_success == False

#     def test_process_blast_output_and_reorient_no_one(self):
#         # Test scenario where the no BLAST hit begins with 1 (start of gene)
#         blast_file = os.path.join(test_data, "SAOMS1_blast_output_no_one.txt")
#         input = os.path.join(test_data, "SAOMS1.fasta")
#         output = os.path.join(test_data, "fake_reoriented.fasta")
#         gene = "terL"
#         blast_success = process_blast_output_and_reorient(
#             input, blast_file, output, gene
#         )
#         assert blast_success == False

#     def test_process_blast_output_and_reorient_correct(self):
#         # Test scenario where the no BLAST hit begins with 1 (start of gene)
#         blast_file = os.path.join(test_data, "SAOMS1_blast_output_correct.txt")
#         input = os.path.join(test_data, "SAOMS1.fasta")
#         output = os.path.join(test_data, "fake_reoriented.fasta")
#         gene = "terL"
#         blast_success = process_blast_output_and_reorient(
#             input, blast_file, output, gene
#         )
#         assert blast_success == True

#     def test_begin_dnaapler(self):
#         # Test begin
#         input = os.path.join(test_data, "SAOMS1.fasta")
#         threads = str(8)
#         gene = "terL"
#         outdir = os.path.join(test_data, "bad_dir")
#         begin_dnaapler(input, outdir, threads, gene)

#     def test_end_dnaapler(self):
#         time = 2324.0
#         end_dnaapler(time)


# class TestEValue(unittest.TestCase):
#     """Tests of Evalue"""

#     def test_evalue_char(self):
#         with self.assertRaises(SystemExit):
#             evalue = "sfsd"
#             check_evalue(evalue)

#     def test_evalue_char_mix(self):
#         with self.assertRaises(SystemExit):
#             evalue = "1e-10t"
#             check_evalue(evalue)

#     def test_evalue_char_int(self):
#         evalue = "5"
#         check_evalue(evalue)

#     def test_evalue_int(self):
#         evalue = 5
#         check_evalue(evalue)

#     def test_evalue_sci(self):
#         evalue = "1e-10"
#         check_evalue(evalue)


# class TestChoiceAutocomplete(unittest.TestCase):
#     """Tests of Choice Autocomplete"""

#     def test_evalue_bad_char(self):
#         # fake values
#         ctx = "1"
#         param = "2"
#         value = "sfsd"
#         with self.assertRaises(click.BadParameter):
#             val = validate_choice_autocomplete(ctx, param, value)

#     def test_evalue_none(self):
#         value = "none"
#         ctx = "1"
#         param = "2"
#         val = validate_choice_autocomplete(ctx, param, value)

#     def test_evalue_mys(self):
#         value = "mystery"
#         ctx = "1"
#         param = "2"
#         val = validate_choice_autocomplete(ctx, param, value)

#     def test_evalue_nearest(self):
#         value = "nearest"
#         ctx = "1"
#         param = "2"
#         val = validate_choice_autocomplete(ctx, param, value)


# # external tools

# """
# taken from tbpore
# """


# class TestExternalTools:
#     @patch.object(
#         ExternalTool,
#         ExternalTool._build_command.__name__,
#         return_value=["mocked", "command", "arg"],
#     )
#     @patch.object(Path, Path.mkdir.__name__)
#     def test___constructor(self, mkdir_mock, build_command_mock):
#         logdir = Path("logs")

#         external_tool = ExternalTool("tool", "input", "output", "params", logdir)

#         assert external_tool.command == ["mocked", "command", "arg"]
#         assert external_tool.command_as_str == "mocked command arg"
#         assert (
#             external_tool.out_log
#             == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.out"
#         )
#         assert (
#             external_tool.err_log
#             == "logs/tool_c238863b32d18040bbf255fa3bf0dc91e9afa268335b56f51abe3c6d1fd83261.err"
#         )

#         build_command_mock.assert_called_once_with("tool", "input", "output", "params")
#         mkdir_mock.assert_called_once_with(parents=True, exist_ok=True)

#     def test___build_command___simple_command(self):
#         expected_escaped_command = ["tool", "param1", "param2", "-o", "out", "-i", "in"]
#         actual_escaped_command = ExternalTool._build_command(
#             "tool", "-i in", "-o out", "param1 param2"
#         )
#         assert expected_escaped_command == actual_escaped_command

#     def test___build_command___single_quote_escaped(self):
#         expected_escaped_command = [
#             "tool",
#             "params",
#             "with",
#             "escaped arg",
#             "-o",
#             "escaped out",
#             "-i",
#             "escaped in",
#         ]
#         actual_escaped_command = ExternalTool._build_command(
#             "tool", "-i 'escaped in'", "-o 'escaped out'", "params with 'escaped arg'"
#         )
#         assert expected_escaped_command == actual_escaped_command

#     def test___build_command___double_quote_escaped(self):
#         expected_escaped_command = [
#             "tool",
#             "params",
#             "with",
#             "escaped arg",
#             "-o",
#             "escaped out",
#             "-i",
#             "escaped in",
#         ]
#         actual_escaped_command = ExternalTool._build_command(
#             "tool", '-i "escaped in"', '-o "escaped out"', 'params with "escaped arg"'
#         )
#         assert expected_escaped_command == actual_escaped_command

#     def test___run(self):
#         logsdir = repo_root.parent.parent / "tests/helpers/logs"
#         logsdir.mkdir(parents=True, exist_ok=True)
#         for file in logsdir.iterdir():
#             file.unlink()

#         python_script = str(repo_root.parent.parent / "tests/helpers/run_test.py")
#         external_tool = ExternalTool(
#             sys.executable,
#             "input",
#             "output",
#             python_script,
#             logsdir,
#         )

#         external_tool.run()

#         out_file = glob.glob(f"{logsdir}/*.out")[0]
#         with open(out_file) as out_file_fh:
#             lines = out_file_fh.readlines()
#             assert lines == ["out\n"]

#         err_file = glob.glob(f"{logsdir}/*.err")[0]
#         with open(err_file) as err_file_fh:
#             lines = err_file_fh.readlines()
#             assert lines == [
#                 "err\n",
#                 f"Command line: {sys.executable} {python_script} output input\n",
#             ]


# class TestFailExternal(unittest.TestCase):
#     """Fail Extenral Tool Test"""

#     def test___run_exit(self):
#         with self.assertRaises(FileNotFoundError):
#             logsdir = repo_root.parent.parent / "tests/helpers/logs"
#             logsdir.mkdir(parents=True, exist_ok=True)
#             for file in logsdir.iterdir():
#                 file.unlink()

#             python_script = str(repo_root.parent.parent / "tests/helpers/run_test.py")
#             external_tool = ExternalTool(
#                 "break_here",
#                 "input",
#                 "output",
#                 python_script,
#                 logsdir,
#             )

#             external_tool.run()
