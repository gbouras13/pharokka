import os
import shutil

import click
import polars as pl
from Bio import SeqIO
from loguru import logger

from .citation import __citation__
from .version import __version__


def get_version():
    return __version__


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
        contig_count = len(list(SeqIO.parse(handle, "fasta")))
    return contig_count


def get_contig_headers(fasta_file) -> pl.Series:
    headers = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.description)
    return pl.Series("header", headers)


# function to touch create a file
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


def rename_file(old_path, new_path):
    """
    Rename a file from old_path to new_path.
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
