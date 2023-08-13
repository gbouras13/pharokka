#!/usr/bin/env python3
import argparse
import os
import sys
from argparse import RawTextHelpFormatter

import input_commands
import plot


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
        help='Output plot file name. ".png" suffix will be added to this automatically. \nWill be output in the pharokka output directory if -o is specified, or in the working directory if --gff andf --genbank are specified/',
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
        help="Number of characters to include in annoation labels before truncation with ellipsis. Must be an integer. Defaults to 20.",
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
        help="Controls the proporition of annotations labelled. Must be a number between 0 and 1 inclusive. \n0 = no annotations, 0.5 = half of the annotations, 1 = all annotations. Defaults to 1. Chosen in order of CDS size.",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_input()

    print("Checking your inputs.")

    try:
        int(args.interval)
    except:
        sys.exit(
            "Error: --interval specified is not an integer. Please check your input and try again. \n"
        )

    try:
        int(args.label_size)
    except:
        sys.exit(
            "Error: --label_size specified is not an integer. Please check your input and try again. \n"
        )

    try:
        int(args.title_size)
    except:
        sys.exit(
            "Error: --title_size specified is not an integer. Please check your input and try again. \n"
        )

    try:
        int(args.dpi)
    except:
        sys.exit(
            "Error: --dpi specified is not an integer. Please check your input and try again. \n"
        )

    try:
        float(args.annotations)
    except:
        sys.exit(
            "Error: --annotations specified is not an float. Please check your input and try again. \n"
        )

    # check if plot already exists
    if args.outdir == "":
        plot_file = str(args.plot_name) + ".png"
    else:
        plot_file = os.path.join(args.outdir, str(args.plot_name) + ".png")

    if args.force == True:
        if os.path.isfile(plot_file) == True:
            os.remove(plot_file)
        else:
            print(
                "\n--force was specified even though the output plot file does not already exist. Continuing \n"
            )
    else:
        if os.path.isfile(plot_file) == True:
            sys.exit(
                "\nOutput plot file already exists and force was not specified. Please specify -f or --force to overwrite the output plot file. \n"
            )

    # flag to see if user provided gff and genbank or output directory
    gff_genbank_flag = True

    if args.gff == "" or args.genbank == "":
        if args.outdir == "":
            sys.exit(
                "\nYou have not specified both a Pharokka gff and a gbk file, or a Pharokka output directory. \nPlease check your input and try again."
            )
        else:
            print("You have specified a Pharokka output directory. \nContinuing.")
            gff_genbank_flag = False

            # check the outdir exists
            if os.path.isdir(args.outdir) == False:
                sys.exit(
                    "\nProvided Pharokka output directory does not exist. Please check your -o or --outdir input. \n"
                )

    else:
        if args.outdir != "":
            print(
                "You have specified Pharokka genbank and gff files and also an output directory.\nThe directory will be ignored. \nContinuing with the gff and genbank file."
            )
        else:
            print(
                "You have specified both Pharokka genbank and gff files. \nContinuing."
            )

    # check the input fasta exists
    input_commands.validate_fasta(args.infile)

    # get num input contigs
    num_contigs = len([1 for line in open(args.infile) if line.startswith(">")])
    if num_contigs > 1:
        print(
            "More than one contig detected in the input file. Only the first phage contig will be plotted. \nContinuing."
        )

    # gff file

    if gff_genbank_flag == True:
        gff_file = args.gff
    else:
        gff_file = os.path.join(args.outdir, args.prefix + ".gff")

    if os.path.isfile(gff_file) == False:
        sys.exit(
            str(gff_file)
            + "was not found. Please check the prefix value `-p` matches the pharokka output directory, \n or use --gff and --genbank to specify the gff and genbank files and try again."
        )

    # Load Genbank file
    gbk_file = os.path.join(args.outdir, args.prefix + ".gbk")

    if gff_genbank_flag == True:
        gbk_file = args.genbank
    else:
        gbk_file = os.path.join(args.outdir, args.prefix + ".gbk")

    # validate gbk exists

    if os.path.isfile(gbk_file) == False:
        sys.exit(
            str(gbk_file)
            + "was not found. Please check the prefix value `-p` matches the pharokka output directory, \n or use --gff and --genbank to specify the gff and genbank files and try again."
        )

    print("All other checked.")
    print("Plotting the phage.")

    plot.create_plot(
        gff_file,
        gbk_file,
        args.interval,
        args.annotations,
        args.title_size,
        args.plot_title,
        args.truncate,
        plot_file,
        args.dpi,
        args.label_size,
        args.label_hypotheticals,
        args.remove_other_features_labels,
    )
