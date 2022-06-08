import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil

v = '1.0.1'


### GLOBAL VARIABLES


def get_input():
	parser = argparse.ArgumentParser(description='phrokka: phage genome annotation piepline', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='where to write the output', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-d', '--database', action="store", help='database directory. If the databases have been install in the default directory, this is not required. Otherwise specify the path',  default='Default')
	parser.add_argument('-f', '--force', help="Overwrites the output directory", action="store_true" )
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, force):
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			shutil.rmtree(output_dir)
		else:
			print("--force was specified even though the outdir does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("Output directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  


	# instantiate outdie

	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	mmseqs_dir = os.path.join(output_dir, "mmseqs/")
	if os.path.isdir(mmseqs_dir) == False:
		os.mkdir(mmseqs_dir)
	return output_dir


def validate_fasta(filename):
	with open(filename, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		if any(fasta):
			print("FASTA checked")
		else:
			sys.exit("Error: Input file is not in the FASTA format.\n")  




