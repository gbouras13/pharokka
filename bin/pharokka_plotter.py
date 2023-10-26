#!/usr/bin/env python3
import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from pathlib import Path

from input_commands import validate_fasta
from loguru import logger
from plot import create_plot
from util import get_version


def get_input():
    """gets input for pharokka_plotter.py
    :return: args
    """
    parser = argparse.ArgumentParser(
        description="pharokka_plotter.py: pharokka plotting function",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        action="store",
        required=True,
        help="Input genome file in FASTA format.",
    )
    parser.add_argument(
        "-n",
        "--plot_name",
        action="store",
        default="pharokka_plot",
        help='Output plot file name. ".png" suffix will be added to this automatically. \nWill be output in the Pharokka output directory if -o is specified, or in the working directory if --gff andf --genbank are specified.',
    )
    parser.add_argument(
        "-o", "--outdir", action="store", default="", help="Pharokka output directory."
    )
    parser.add_argument("--gff", action="store", default="", help="Pharokka gff.")
    parser.add_argument(
        "--genbank", action="store", default="", help="Pharokka genbank."
    )
    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        help="Prefix used to create pharokka output. Will default to pharokka.",
        default="pharokka",
    )
    parser.add_argument(
        "-t", "--plot_title", action="store", default="Phage", help="Plot name."
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
    logger.info("Running pharokka_plotter.py to plot your phage.")
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
            f"--annotations {args.annotations} specified is not an float. Please check your input and try again."
        )

    # check if plot already exists
    if args.outdir == "":
        plot_file = str(args.plot_name) + ".png"
        svg_plot_file = str(args.plot_name) + ".svg"
    else:
        plot_file = os.path.join(args.outdir, f"{args.plot_name}.png")
        svg_plot_file = os.path.join(args.outdir, f"{args.plot_name}.svg")

    if args.force == True:
        if os.path.isfile(plot_file) == True:
            os.remove(plot_file)
        else:
            logger.warning(
                "--force was specified even though the output plot file does not already exist."
            )
            logger.warning("Continuing")
    else:
        if os.path.isfile(plot_file) == True:
            logger.error(
                f"Output plot file {plot_file} already exists and force was not specified. Please specify -f or --force to overwrite the output plot file."
            )

    if args.force == True:
        if os.path.isfile(svg_plot_file) == True:
            os.remove(svg_plot_file)
        else:
            logger.warning(
                "--force was specified even though the output plot file does not already exist."
            )
            logger.warning("Continuing")
    else:
        if os.path.isfile(svg_plot_file) == True:
            logger.error(
                f"Output plot file {svg_plot_file} already exists and force was not specified. Please specify -f or --force to overwrite the output plot file."
            )

    # flag to see if user provided gff and genbank or output directory
    gff_genbank_flag = True

    if args.gff == "" or args.genbank == "":
        if args.outdir == "":
            logger.error(
                "You have not specified both a Pharokka gff and a gbk file, or a Pharokka output directory. Please check your input and try again."
            )
        else:
            logger.info("You have specified a Pharokka output directory.")
            gff_genbank_flag = False

            # check the outdir exists
            if os.path.isdir(args.outdir) == False:
                logger.error(
                    f"Provided Pharokka output directory {args.outdir} does not exist. Please check -o/--outdir."
                )

    else:
        if args.outdir != "":
            logger.info(
                "You have specified Pharokka genbank and gff files and also an output directory."
            )
            logger.info(
                "The directory will be ignored. Continuing with the gff and genbank file."
            )
        else:
            logger.info("You have specified both Pharokka genbank and gff files.")

    # check the input fasta exists
    validate_fasta(args.infile)

    # get num input contigs
    num_contigs = len([1 for line in open(args.infile) if line.startswith(">")])
    if num_contigs > 1:
        logger.warning(
            "More than one contig detected in the input file. Only the first phage contig will be plotted."
        )

    # gff file

    if gff_genbank_flag == True:
        gff_file = args.gff
    else:
        gff_file = os.path.join(args.outdir, args.prefix + ".gff")

    if os.path.isfile(gff_file) == False:
        logger.error(
            f" {gff_file} was not found. Please check the prefix value `-p` matches the pharokka output directory, \n or use --gff and --genbank to specify the gff and genbank files and try again."
        )

    # Load Genbank file
    gbk_file = os.path.join(args.outdir, args.prefix + ".gbk")

    if gff_genbank_flag == True:
        gbk_file = args.genbank
    else:
        gbk_file = os.path.join(args.outdir, args.prefix + ".gbk")

    # validate gbk exists

    if os.path.isfile(gbk_file) == False:
        logger.error(
            f" {gbk_file} was not found. Please check the prefix value `-p` matches the pharokka output directory, \n or use --gff and --genbank to specify the gff and genbank files and try again."
        )

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
    logger.info("Plotting the phage.")

    create_plot(
        gff_file,
        gbk_file,
        args.interval,
        args.annotations,
        args.title_size,
        args.plot_title,
        args.truncate,
        plot_file,
        svg_plot_file,
        args.dpi,
        args.label_size,
        args.label_hypotheticals,
        args.remove_other_features_labels,
        label_force_list,
    )
