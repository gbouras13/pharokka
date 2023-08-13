#!/usr/bin/env python3
import argparse
import os
from argparse import RawTextHelpFormatter

import databases


def get_db_input():
    parser = argparse.ArgumentParser(
        description="script to download required phrog databases",
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
        help="Database Directory - will be created and must be specificed if -d is not used. ",
    )
    args = parser.parse_args()
    return args


args = get_db_input()

if args.default == True:
    db_dir = os.path.join(os.path.dirname(__file__), "../", "databases/")
    databases.instantiate_install(db_dir)
else:
    if args.outdir == None:
        print("--outdir was not specified.")
        print("Downloading databases to the default directory anyway.")
        db_dir = os.path.join(os.path.dirname(__file__), "../", "databases/")
    else:
        db_dir = args.outdir
    databases.instantiate_install(db_dir)
