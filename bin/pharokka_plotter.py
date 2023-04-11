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
    parser.add_argument('-o', '--output', action="store", required = True, help='Output directory.')
    parser.add_argument('-p', '--prefix', action="store", help='Prefix used to create pharokka output. Will default to pharokka.',  default='pharokka')
    parser.add_argument('--plot_name', action="store",  default='Phage', help='Plot name.')
    parser.add_argument('--interval', action="store",  default='5000', help='axis ticks interval.')
    parser.add_argument('--title_size', action="store",  default='20', help='Controls title size.')
    parser.add_argument('--sparse', action="store",  default='1', help='Plot fewer annotations if plot is crowded. 2 = half, 3 = third etc. Must be an integer')
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = get_input()
    print("Plotting the phage.")

    # get num input contigs
    num_contigs = len([1 for line in open(args.infile) if line.startswith(">")])
    if num_contigs > 1:
        print("More than one contig detected in the input file. Only the first phage contig will be plotted. \nContinuing.")
    plot.create_plot(args.output, args.prefix, args.plot_name, args.interval, args.sparse, args.title_size)

