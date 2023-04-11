#!/usr/bin/env python3
import argparse
import os
from argparse import RawTextHelpFormatter
import plot



def get_input():
    """gets input for pharokka_plotter.py
    :return: args
    """
    parser = argparse.ArgumentParser(description='pharokka_plotter.py: pharokka plotting function', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", required = True, help='Input genome file in FASTA format.')
    parser.add_argument('-d', '--directory', action="store", required = True, help='Pharokka output directory.')
    parser.add_argument('-o', '--outfile', action="store", default='pharokka_plot.png', help='Output png file name')
    parser.add_argument('-p', '--prefix', action="store", help='Prefix used to create pharokka output. Will default to pharokka.',  default='pharokka')
    parser.add_argument('-t', '--plot_title', action="store",  default='Phage', help='Plot name.')
    parser.add_argument('--title_size', action="store",  default='20', help='Controls title size.')
    parser.add_argument('--label_size', action="store",  default='8', help='Controls annotation label size.')
    parser.add_argument('--interval', action="store",  default='5000', help='Axis tick interval.')
    parser.add_argument('--truncate', action="store",  default='20', help='Number of characters to include in annoation labels before truncation with ellipsis.')
    parser.add_argument('--dpi', action="store",  default='600', help='Resultion (dots per inch). Defaults to 600.')
    parser.add_argument('--annotations', action="store",  default='1', help='Controls the proporition of annotations labelled (between 0 and 1). 0 = no annotations, 1 = all annotations. Plots in order of CDS size.')
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = get_input()
    print("Plotting the phage.")

    # get num input contigs
    num_contigs = len([1 for line in open(args.infile) if line.startswith(">")])
    if num_contigs > 1:
        print("More than one contig detected in the input file. Only the first phage contig will be plotted. \nContinuing.")
    plot.create_plot(args.directory, args.prefix, args.interval, args.annotations, args.title_size, args.plot_title, args.truncate, args.outfile, args.dpi, args.label_size)

