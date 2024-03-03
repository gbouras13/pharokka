#!/usr/bin/env python3
import argparse
import os
import shutil
import sys
from argparse import RawTextHelpFormatter
from pathlib import Path

from loguru import logger
from plot import create_single_plot
from pycirclize.parser import Genbank
from util import get_version


def get_input():
    """gets input for pharokka_plotter.py
    :return: args
    """
    parser = argparse.ArgumentParser(
        description="pharokka_multiplotter.py: pharokka plotting function for muliple phages",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-g",
        "--genbank",
        action="store",
        required=True,
        help="Input genbank file from Pharokka.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Pharokka output directory.",
        required=True,
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output plot file.", action="store_true"
    )
    parser.add_argument(
        "--label_hypotheticals",
        help="Flag to label  hypothetical or unknown proteins. By default these are not labelled.",
        action="store_true",
    )
    parser.add_argument(
        "--remove_other_features_labels",
        help="Flag to remove labels for tRNA/tmRNA/CRISPRs. By default these are labelled. \nThey will still be plotted in black.",
        action="store_true",
    )
    parser.add_argument(
        "--title_size",
        action="store",
        default="20",
        help="Controls title size. Must be an integer. Defaults to 20.",
    )
    parser.add_argument(
        "--label_size",
        action="store",
        default="8",
        help="Controls annotation label size. Must be an integer. Defaults to 8.",
    )
    parser.add_argument(
        "--interval",
        action="store",
        default="5000",
        help="Axis tick interval. Must be an integer. Must be an integer. Defaults to 5000.",
    )
    parser.add_argument(
        "--truncate",
        action="store",
        default="20",
        help="Number of characters to include in annoation labels before truncation with ellipsis. \nMust be an integer. Defaults to 20.",
    )
    parser.add_argument(
        "--dpi",
        action="store",
        default="600",
        help="Resultion (dots per inch). Must be an integer. Defaults to 600.",
    )
    parser.add_argument(
        "--annotations",
        action="store",
        default="1",
        help="Controls the proporition of annotations labelled. Must be a number between 0 and 1 inclusive. \n0 = no annotations, 0.5 = half of the annotations, 1 = all annotations. \nDefaults to 1. Chosen in order of CDS size.",
    )
    parser.add_argument(
        "-t", "--plot_title", action="store", default="Phage", help="Plot name."
    )
    parser.add_argument(
        "--label_ids",
        action="store",
        default="",
        help="Text file with list of CDS IDs (from gff file) that are guaranteed to be labelled.",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_input()
    # preamble
    logger.add(lambda _: sys.exit(1), level="ERROR")
    logger.info(f"Starting Pharokka v{get_version()}")
    logger.info("Running pharokka_multiplotter.py to plot your phages.")
    logger.info("Command executed: {}", args)
    logger.info("Repository homepage is https://github.com/gbouras13/pharokka")
    logger.info("Written by George Bouras: george.bouras@adelaide.edu.au")
    logger.info("Checking your inputs.")

    try:
        int(args.interval)
    except:
        logger.error(
            f"--interval {args.interval} specified is not an integer. Please check your input and try again."
        )

    try:
        int(args.label_size)
    except:
        logger.error(
            f"--label_size {args.label_size} specified is not an integer. Please check your input and try again."
        )

    try:
        int(args.title_size)
    except:
        logger.error(
            f"--title_size {args.title_size} specified is not an integer. Please check your input and try again."
        )

    try:
        int(args.dpi)
    except:
        logger.error(
            f"--dpi {args.dpi} specified is not an integer. Please check your input and try again."
        )

    try:
        float(args.annotations)
    except:
        logger.error(
            f"--annotations {args.annotations} specified is not a float. Please check your input and try again."
        )

    if args.force == True:
        if os.path.exists(args.outdir) == True:
            logger.info(f"Removing {args.outdir} as --force was specified.")
            shutil.rmtree(args.outdir)
        else:
            logger.warning(
                f"--force was specified even though the output plot directory {args.outdir} does not already exist."
            )
            logger.warning("Continuing")
    else:
        if os.path.exists(args.outdir) == True:
            logger.error(
                f"Output directoey {args.outdir} already exists and force was not specified. Please specify -f or --force to overwrite the output directory."
            )

    # instantiate outdir
    if Path(args.outdir).exists() is False:
        Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # check label_ids

    # list of all IDs that need to be labelled from file
    label_force_list = []

    if args.label_ids != "":
        logger.info(
            f"You have specified a file {args.label_ids} containing a list of CDS IDs to force label."
        )
        # check if it is a file
        if os.path.isfile(args.label_ids) == False:
            logger.error(f"{args.label_ids} was not found.")
        # check if it contains text
        try:
            # Open the file in read mode
            with open(Path(args.label_ids), "r") as file:
                # Read the first character
                first_char = file.read(1)

                # read all the labels
                with open(Path(args.label_ids)) as f:
                    ignore_dict = {x.rstrip().split()[0] for x in f}
                # label force list
                label_force_list = list(ignore_dict)

        except FileNotFoundError:
            logger.warning(
                f"{args.label_id} contains no text. No contigs will be ignored"
            )

    logger.info("All files checked.")

    # single threaded plots
    threads = 1

    gbk = Genbank(args.genbank)

    # gets all contigs and seqs
    gb_seq_dict = gbk.get_seqid2seq()

    gb_size_dict = gbk.get_seqid2size()

    contig_count = len(gb_seq_dict)

    # gets all features - will get all regardless of type (tRNA etc from pharokka)
    gb_feature_dict = gbk.get_seqid2features()

    # check label_ids
    # list of all IDs that need to be labelled from file
    label_force_list = []

    label_ids = args.label_ids

    if label_ids != "":
        logger.info(
            f"You have specified a file {label_ids} containing a list of CDS IDs to force label."
        )
        # check if it is a file
        if Path(label_ids).exists() is False:
            logger.error(f"{label_ids} was not found.")
        # check if it contains text
        try:
            # Open the file in read mode
            with open(Path(label_ids), "r") as file:
                # Read the first character
                # will error if file is empty
                first_char = file.read(1)

                # read all the labels
                with open(Path(label_ids)) as f:
                    ignore_dict = {x.rstrip().split()[0] for x in f}
                # label force list
                label_force_list = list(ignore_dict)

        except FileNotFoundError:
            logger.warning(f"{label_ids} contains no text. No contigs will be ignored")

    # if there is 1 contig, then all the parameters will apply

    for contig_id, contig_sequence in gb_seq_dict.items():
        logger.info(f"Plotting {contig_id}")

        create_single_plot(
            contig_id,
            contig_sequence,
            contig_count,
            gb_size_dict,
            gb_feature_dict,
            gbk,
            int(args.interval),
            float(args.annotations),
            int(args.title_size),
            args.plot_title,
            int(args.truncate),
            args.outdir,
            int(args.dpi),
            float(args.label_size),
            args.label_hypotheticals,
            args.remove_other_features_labels,
            label_force_list,
        )
