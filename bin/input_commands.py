import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil
from version import __version__

v = __version__


### GLOBAL VARIABLES


def get_input():
	parser = argparse.ArgumentParser(description='pharokka: fast phage annotation program', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='Input genome file in fasta format.',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Directory to write the output to.', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-d', '--database', action="store", help='Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.',  default='Default')
	parser.add_argument('-t', '--threads', help="Number of threads for mmseqs and hhsuite. Defaults to 1.", action="store", default = str(1))
	parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
	parser.add_argument('-p', '--prefix', action="store", help='Prefix for output files. This is not required.',  default='Default')
	parser.add_argument('-l', '--locustag', action="store", help='User specified locus tag for the gff/gbk files. This is not required. A random locus tag will be generated instead.',  default='Default')
	parser.add_argument('-g', '--gene_predictor', action="store", help='User specified gene predictor. Use "-g phanotate" or "-g prodigal". Defaults to phanotate (not required unless prodigal is desired).',  default='phanotate' )
	parser.add_argument('-m', '--meta', help='meta mode for metavirome input samples', action="store_true")
	parser.add_argument('-c', '--coding_table', help='translation table for prodigal. Defaults to 11. Experimental only.', action="store", default = "11")
	parser.add_argument('-e', '--evalue', help='E-value threshold for mmseqs2 PHROGs database search. Defaults to 1E-05.', action="store", default = "1E-05")
	parser.add_argument('-te', '--terminase', help='Runs terminase large subunit re-orientation mode.', action="store_true")
	parser.add_argument('-s', '--strand', help='Strand of terminase large subunit. Must be "pos" or "neg".', action="store", default = "nothing")
	parser.add_argument('-ts', '--terminase_start', help='Start coordinate of the terminase large subunit.', action="store", default = "0")
	parser.add_argument('-V', '--version', help='Version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir, meta, gene_predictor, force):
	# remove outdir on force
	if force == True:
		if os.path.isdir(output_dir) == True:
			shutil.rmtree(output_dir)
		else:
			print("\n--force was specified even though the outdir does not already exist. Continuing \n")
	else:
		if os.path.isdir(output_dir) == True:
			sys.exit("\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  


	# instantiate outdir
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	mmseqs_dir = os.path.join(output_dir, "mmseqs/")
	if os.path.isdir(mmseqs_dir) == False:
		os.mkdir(mmseqs_dir)
	vfdb_dir = os.path.join(output_dir, "vfdb/")
	if os.path.isdir(vfdb_dir) == False:
		os.mkdir(vfdb_dir)
	CARD_dir = os.path.join(output_dir, "CARD/")
	if os.path.isdir(CARD_dir) == False:
		os.mkdir(CARD_dir)

	# tmp dir for  meta mode trnascan and phanotate
	input_tmp_dir = os.path.join(output_dir, "input_split_tmp/")
	if meta == True:
		if os.path.isdir(input_tmp_dir) == False:
			os.mkdir(input_tmp_dir)

	return output_dir

def validate_fasta(filename):
	with open(filename, "r") as handle:
		fasta = SeqIO.parse(handle, "fasta")
		print("Checking Input FASTA")
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

def validata_meta(filepath_in, meta):
	num_fastas = len([1 for line in open(filepath_in) if line.startswith(">")])
	if meta == True:
		if num_fastas < 2:
			sys.exit("Error: -m meta mode specified when the input file only contains 1 contig. Please re-run without specifying -m. \n")  
	else:
		if num_fastas > 1:
			print("More than one contig detected in the input file. Re-running pharokka with -m meta mode is recommended. \n Continuing.")

def validate_strand(strand):
	if strand != "pos" and strand != "neg":
		sys.exit("Error: terminase strand was incorrectly specified. It should be either 'pos' or 'neg'. Please check your input and try again. \n")  

def validate_terminase_start(terminase_start):
    try:
        int(terminase_start)
    except:
        sys.exit("Error: terminase start coordinate specified is not an integer. Please check your input and try again. \n") 


def validate_terminase(filepath_in, strand, terminase_start):
	if strand == "nothing":
		sys.exit("Error: you specified -te to reorient your phage to begin with the terminase large subunit, but didn't specify its strand with -s. Please check your input and try again. \n") 
	if terminase_start == "0":
		sys.exit("Error: you specified -te to reorient your phage to begin with the terminase large subunit, but didn't specify its start coordinate with -ts. Please check your input and try again. \n") 
	validate_strand(strand)
	validate_terminase_start(terminase_start)
	num_fastas = len([1 for line in open(filepath_in) if line.startswith(">")])
	if num_fastas > 1:
		sys.exit("Error: To reorient your phage genome to begin with the terminase large subunit, you can only input 1 phage genome. Multiple contigs were detected. Please try again. \n")  
