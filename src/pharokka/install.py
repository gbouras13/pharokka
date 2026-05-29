"""
Entry point for `pharokka install` – download/install databases.

Adapted from bin/install_databases.py with relative package imports.
"""

import argparse
import os
import sys
from argparse import RawTextHelpFormatter

from loguru import logger

from .databases import instantiate_install


def get_db_input():
    parser = argparse.ArgumentParser(
        description="script to download required Pharokka databases",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--default",
        action="store_true",
        help="Determines whether you want databases stored in the default location.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Database Directory. It will be created and must be specified if -d is not used. ",
    )
    args = parser.parse_args()
    return args


def _default_db_dir():
    """Return the databases/ directory at the project root (3 levels above this file)."""
    _here = os.path.dirname(os.path.realpath(__file__))  # src/pharokka/
    _src = os.path.dirname(_here)                         # src/
    _root = os.path.dirname(_src)                         # project root
    return os.path.join(_root, "databases/")


def main():
    logger.add(lambda _: sys.exit(1), level="ERROR")
    args = get_db_input()

    if args.default is True:
        db_dir = _default_db_dir()
        instantiate_install(db_dir)
    else:
        if args.outdir is None:
            logger.info("--outdir was not specified.")
            logger.info("Downloading databases to the default directory anyway.")
            db_dir = _default_db_dir()
        else:
            db_dir = args.outdir
        instantiate_install(db_dir)


if __name__ == "__main__":
    main()
