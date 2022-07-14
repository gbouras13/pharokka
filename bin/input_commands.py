import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil

v = '0.1.6'


### GLOBAL VARIABLES


def get_input():
	parser = argparse.ArgumentParser(description='pharokka: phage genome annotation piepline', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='Input file in fasta format.',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to.', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-d', '--database', action="store", help='Database directory. If the databases have been install in the default directory, this is not required. Otherwise specify the path.',  default='Default')
	parser.add_argument('-t', '--threads', help="Number of threads for mmseqs and hhsuite. Defaults to 1.", action="store", default = str(1))
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required',  default='Default')
	parser.add_argument('-l', '--locustag', action="store", help='User specified locus tag for the gff/gbk files. This is not required. A random locus tag will be generated instead.',  default='Default')
	parser.add_argument('-g', '--gene_predictor', action="store", help='User specified gene predictor. Use " -g phanotate" or "-g prodigal". Defaults to phanotate.',  default='phanotate' )
	parser.add_argument('-m', '--meta', help='Metagenomic option for Prodigal', action="store_true")
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, force):
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			shutil.rmtree(output_dir)
		else:
			print("\n--force was specified even though the outdir does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  


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

def validate_gene_predictor(gene_predictor):
	if gene_predictor == "phanotate":
		print("Phanotate will be used for gene prediction")
	elif gene_predictor == "prodigal":
		print("Prodigal will be used for gene prediction")
	else:
		sys.exit("Error: gene predictor was incorrectly specified. Please use 'phanotate' or 'prodigal' .\n")  



