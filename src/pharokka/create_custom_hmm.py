#!/usr/bin/env python3

"""
Script to create hmm profile with pyhmmer

Note: all MSAs must be in FASTA format inside the directory provided with -i and labelled with only 1 full stop e.g. "name.msa"
"""

import argparse
import os
import shutil
from argparse import RawTextHelpFormatter
from pathlib import Path

import pyhmmer
from loguru import logger

from .util import get_version

alphabet = pyhmmer.easel.Alphabet.amino()
background = pyhmmer.plan7.Background(alphabet)


def get_input():
    """gets input for create_custom_hmm
    :return: args
    """
    parser = argparse.ArgumentParser(
        description="pharokka create-hmm: Creates HMMs from FASTA formatted MSAs with PyHMMER for use with Pharokka v1.4.0 and higher.",
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

    # The logger.error → sys.exit(1) sink is registered once per process
    # by the entry-point dispatcher (cli.main or pharokka_scripts._legacy_shim).
    logger.info(f"Starting pharokka v{get_version()} - create_hmms")

    MSA_dir = args.indir
    HMM_dir = args.outdir

    if args.force:
        if os.path.isdir(HMM_dir):
            logger.info(
                f"Removing output directory {HMM_dir} as -f or --force was specified."
            )
            shutil.rmtree(HMM_dir)
        elif os.path.isfile(HMM_dir):
            logger.info(
                f"Removing output file {HMM_dir} as -f or --force was specified."
            )
            os.remove(HMM_dir)
        else:
            logger.info(
                f"--force was specified even though the output directory {HMM_dir} does not already exist. Continuing."
            )
    else:
        if os.path.isdir(HMM_dir) or os.path.isfile(HMM_dir):
            logger.error(
                f"The output directory {HMM_dir} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )

    if not os.path.exists(HMM_dir):
        os.mkdir(HMM_dir)

    logger.info(
        f"Creating HMMs in the directory {HMM_dir} from MSAs in the directory {MSA_dir}."
    )

    file_list = os.listdir(MSA_dir)

    for file in file_list:
        if file.startswith("."):
            continue
        else:
            if is_fasta_msa(f"{MSA_dir}/{file}"):
                with pyhmmer.easel.MSAFile(
                    f"{MSA_dir}/{file}", digital=True, alphabet=alphabet
                ) as msa_file:
                    msa = msa_file.read()
                root, _ = os.path.splitext(file)
                name = root
                msa.name = name
                builder = pyhmmer.plan7.Builder(alphabet)
                background = pyhmmer.plan7.Background(alphabet)
                hmm, _, _ = builder.build_msa(msa, background)
                with open(f"{HMM_dir}/{name}.hmm", "wb") as output_file:
                    hmm.write(output_file)
            else:
                logger.warning(
                    f"{MSA_dir}/{file} does not seem to be a FASTA formatted MSA. Skipping."
                )

    hmms = []
    HMM_dir = Path(HMM_dir)
    hmm_file_list = os.listdir(HMM_dir)

    for file_name in hmm_file_list:
        f = f"{HMM_dir}/{file_name}"
        with pyhmmer.plan7.HMMFile(f) as hmm_file:
            hmm = hmm_file.read()
        hmms.append(hmm)

    pyhmmer.hmmer.hmmpress(hmms, f"{HMM_dir}/{args.prefix}")

    logger.info("HMM creation complete.")
    logger.info(
        f"The combined file you will need to run with pharokka run --custom_hmm is {HMM_dir}/{args.prefix}.h3m"
    )


def is_fasta_msa(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    return len([line for line in lines if line.startswith(">")]) > 1


if __name__ == "__main__":
    main()
