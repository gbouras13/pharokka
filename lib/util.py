import os
import shutil

import click
from loguru import logger

from lib.citation import __citation__
from lib.version import __version__


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
