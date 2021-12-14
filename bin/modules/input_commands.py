import argparse
import os
import VERSION
from argparse import RawTextHelpFormatter
import datetime

v = VERSION.get_version()


### GLOBAL VARIABLES


def get_input():
	usage = 'phrokka ...'
	parser = argparse.ArgumentParser(description='phrokka: phage genome annotation piepline', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='where to write the output', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir):
	# if the output directory already exists add the date and time on the end to make unique dir
	if os.path.isdir(output_dir) == True:
		output_dir = output_dir + str(datetime.datetime.now().strftime('_%Y%m%d_%H%M_%S%f'))
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	mmseqs_dir = os.path.join(output_dir, "mmseqs/")
	if os.path.isdir(mmseqs_dir) == False:
		os.mkdir(mmseqs_dir)





