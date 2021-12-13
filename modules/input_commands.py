import argparse
import os
import sys
from argparse import RawTextHelpFormatter

ROOT_DIR = os.path.dirname(os.path.abspath("VERSION"))
version_file = open(os.path.join(ROOT_DIR, 'VERSION'))
v = version_file.read().strip()

def get_input():
	usage = 'phrokka ...'
	parser = argparse.ArgumentParser(description='phrokka: phage genome annotation piepline', formatter_class=RawTextHelpFormatter)
	parser.add_argument( '-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='where to write the output [stdout]', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir):
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	mmseqs_dir = os.path.join(output_dir, "mmseqs/")
	if os.path.isdir(mmseqs_dir) == False:
		os.mkdir(mmseqs_dir)





