#!/usr/bin/env python3

"""
Script to create hmm profile with pyhmmer

create_custom_hmm.py -i <directory of MSAs> -o <directory with HMMs>

Note: all MSAs must be in FASTA format inside the directory provided with -i and labelled with only 1 full stop e.g. "name.msa"

"""

import os
import shutil
from pathlib import Path

import pyhmmer
from util import get_version

alphabet = pyhmmer.easel.Alphabet.amino()
background = pyhmmer.plan7.Background(alphabet)
import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from pathlib import Path

from loguru import logger


def get_input():
    """gets input for create_custom_hmm.py
    :return: args
    """
    parser = argparse.ArgumentParser(
        description="create_custom_hmm.py: Creates HMMs from FASTA formatted MSAs with PyHMMER for use with Pharokka v1.4.0 and higher.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--indir",
        action="store",
        help="Input directory containing FASTA formatted MSAs (Multiple Sequence Alignments).",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        default="",
        help="Output directory to store HMM profiles.",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        help="Prefix used to name the combined HMM file. The relevant file will be prefix.h3m",
        default="custom_db",
    )

    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    parser.add_argument(
        "-V",
        "--version",
        help="Print pharokka Version",
        action="version",
        version=get_version(),
    )
    args = parser.parse_args()

    return args


def main():
    args = get_input()

    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Starting pharokka v{get_version()} - create_hmms.py")

    MSA_dir = args.indir
    HMM_dir = args.outdir

    #### force
    if args.force == True:
        if os.path.isdir(HMM_dir) == True:
            logger.info(
                f"Removing output directory {HMM_dir} as -f or --force was specified."
            )
            shutil.rmtree(HMM_dir)
        elif os.path.isfile(HMM_dir) == True:
            logger.info(
                f"Removing output file {HMM_dir} as -f or --force was specified."
            )
            os.remove(HMM_dir)
        else:
            logger.info(
                f"--force was specified even though the output directory {HMM_dir} does not already exist. Continuing."
            )
    else:
        if os.path.isdir(HMM_dir) == True or os.path.isfile(HMM_dir) == True:
            logger.error(
                f"The output directory {HMM_dir} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )

    # Check if the directory already exists and make the dir
    if not os.path.exists(HMM_dir):
        # Create the output directory
        os.mkdir(HMM_dir)

    logger.info(
        f"Creating HMMs in the directory {HMM_dir} from MSAs in the directory {MSA_dir}."
    )

    # Get a list of all files in the directory
    file_list = os.listdir(MSA_dir)

    # loop over each PHROG
    for file in file_list:
        # check if MSA
        if is_fasta_msa(f"{MSA_dir}/{file}"):
            # read in each msa
            with pyhmmer.easel.MSAFile(
                f"{MSA_dir}/{file}", digital=True, alphabet=alphabet
            ) as msa_file:
                msa = msa_file.read()
            # split the file into root and suffix
            root, _ = os.path.splitext(file)
            name = root
            # convert to bytes
            msa.name = name.encode("utf-8")
            # build the HMM
            builder = pyhmmer.plan7.Builder(alphabet)
            background = pyhmmer.plan7.Background(alphabet)
            hmm, _, _ = builder.build_msa(msa, background)
            with open(f"{HMM_dir}/{name}.hmm", "wb") as output_file:
                hmm.write(output_file)
        else:
            logger.warning(
                f"{MSA_dir}/{file} does not seem to be a FASTA formatted MSA. Skipping."
            )

    # to concatenate all hmms

    hmms = []

    # Specify the directory path
    HMM_dir = Path(HMM_dir)

    # Get a list of all files in the directory
    hmm_file_list = os.listdir(HMM_dir)

    # reads and saves the hmms
    for file_name in hmm_file_list:
        f = f"{HMM_dir}/{file_name}"
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()
        hmms.append(hmm)

    # writes all out together to .h3m, .h3p, .h3i, .h3f files prefixed "prefix"
    pyhmmer.hmmer.hmmpress(hmms, f"{HMM_dir}/{args.prefix}")

    logger.info(f"HMM creation complete.")
    logger.info(
        f"The combined file you will need to run with pharokka.py --custom_hmm is {HMM_dir}/{args.prefix}.h3m"
    )


def is_fasta_msa(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    return len([line for line in lines if line.startswith(">")]) > 1


if __name__ == "__main__":
    main()
