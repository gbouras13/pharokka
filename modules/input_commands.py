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
	parser.add_argument('-o', '--outdir', action="store", help='where to write the output [stdout]', default=sys.stdout)
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_outdir():
	out_dir = "output/"
	if os.path.isdir(out_dir) == False:
		os.mkdir(out_dir)

