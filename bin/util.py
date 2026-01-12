import os
import shutil

import click
import pandas as pd
from Bio import SeqIO
from citation import __citation__
from loguru import logger
from version import __version__


def get_version():
    version = __version__
    return version


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        with open(log, "a") as lo:
            lo.write(msg)


def print_citation():
    echo_click(__citation__)


log_fmt = (
    "[<green>{time:YYYY-MM-DD HH:mm:ss}</green>] <level>{level: <8}</level> | "
    "<level>{message}</level>"
)


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)


def count_contigs(input_fasta) -> int:
    """
    counts the number of contigs in the FASTA input file
    return: contig_count, int the number of contigs
    """

    with open(input_fasta, "r") as handle:
        # Check the number of records
        contig_count = len(list(SeqIO.parse(handle, "fasta")))

    return contig_count


def get_contig_headers(fasta_file) -> pd.Series:
    headers = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.description)
    headers_series = pd.Series(headers)
    return headers_series


# function to touch create a file
# https://stackoverflow.com/questions/12654772/create-empty-file-using-python
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


def rename_file(old_path, new_path):
    """
    Rename a file from old_path to new_path.

    Parameters:
    - old_path (str): The current path to the file.
    - new_path (str): The new desired path and filename.

    Raises:
    - FileNotFoundError: If the file at old_path does not exist.
    - FileExistsError: If a file already exists at new_path.
    - OSError: For other OS-related errors.
    """
    if not os.path.exists(old_path):
        logger.error(f"File '{old_path}' does not exist.")
    if os.path.exists(new_path):
        logger.error(f"A file already exists at '{new_path}'.")

    os.rename(old_path, new_path)


# Split attributes/description column from pyrodigal gff into key-value pairs
def parse_attributes(attr_str):
    if not attr_str:
        return {}
    pairs = attr_str.split(";")
    return dict(pair.split("=", 1) for pair in pairs if "=" in pair)
