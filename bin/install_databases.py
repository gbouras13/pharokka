#!/usr/bin/env python3
import argparse
import os
import sys
from argparse import RawTextHelpFormatter

from databases import instantiate_install
from loguru import logger


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


logger.add(lambda _: sys.exit(1), level="ERROR")
args = get_db_input()

if args.default == True:
    db_dir = os.path.join(os.path.dirname(__file__), "../", "databases/")
    instantiate_install(db_dir)
else:
    if args.outdir == None:
        logger.info("--outdir was not specified.")
        logger.info("Downloading databases to the default directory anyway.")
        db_dir = os.path.join(os.path.dirname(__file__), "../", "databases/")
    else:
        db_dir = args.outdir

    instantiate_install(db_dir)
