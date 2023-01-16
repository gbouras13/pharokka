import argparse
import os
import sys
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import shutil
from version import __version__
import subprocess as sp


v = __version__


### GLOBAL VARIABLES


def get_input():
	parser = argparse.ArgumentParser(description='pharokka: fast phage annotation program', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='Input genome file in fasta format.')
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
	parser.add_argument('--terminase', help='Runs terminase large subunit re-orientation mode. Single genome input only and requires -s and -te to be specified.', action="store_true")
	parser.add_argument('--terminase_strand', help='Strand of terminase large subunit. Must be "pos" or "neg".', action="store", default = "nothing")
	parser.add_argument('--terminase_start', help='Start coordinate of the terminase large subunit.', action="store", default = "nothing")
	parser.add_argument('-V', '--version', help='Print pharokka Version', action='version', version=v)
	parser.add_argument('--citation', help='Print pharokka Citation', action="store_true")
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
			sys.exit("ERROR: -m meta mode specified when the input file only contains 1 contig. Please re-run without specifying -m. \n")  
	else:
		if num_fastas > 1:
			print("More than one contig detected in the input file. Re-running pharokka with -m meta mode is recommended unless this is a fragmented isolate genome. \nContinuing.")

def validate_strand(strand):
	if strand != "pos" and strand != "neg":
		sys.exit("Error: terminase strand was incorrectly specified. It should be either 'pos' or 'neg'. Please check your input and try again. \n")  

def validate_terminase_start(terminase_start):
    try:
        int(terminase_start)
    except:
        sys.exit("Error: terminase start coordinate specified is not an integer. Please check your input and try again. \n") 


def validate_terminase(filepath_in, terminase_strand, terminase_start):
	if terminase_strand == "nothing":
		sys.exit("Error: you specified -te to reorient your phage to begin with the terminase large subunit, but didn't specify its strand with -s. Please check your input and try again. \n") 
	if terminase_start == "nothing":
		sys.exit("Error: you specified -te to reorient your phage to begin with the terminase large subunit, but didn't specify its start coordinate with -ts. Please check your input and try again. \n") 
	validate_strand(terminase_strand)
	validate_terminase_start(terminase_start)
	num_fastas = len([1 for line in open(filepath_in) if line.startswith(">")])
	if num_fastas > 1:
		sys.exit("Error: To reorient your phage genome to begin with the terminase large subunit, you can only input 1 phage genome contig. Multiple contigs were detected. Please try again. \n")  


#######
# dependencies
#######


def check_dependencies(logger):
	"""Checks the dependencies and versions
    :return:
    """
	#############
	# phanotate
	#############
	try:
		process = sp.Popen(["phanotate.py", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	except:
		sys.exit("Phanotate not found. Please reinstall pharokka.")
	phan_out, _ = process.communicate()
	phanotate_out = phan_out.decode().strip()
	phanotate_major_version = int(phanotate_out.split('.')[0])
	phanotate_minor_version = int(phanotate_out.split('.')[1])
	phanotate_minorest_version = phanotate_out.split('.')[2]

	print("Phanotate version found is v" + str(phanotate_major_version) +"." + str(phanotate_minor_version) +"."+phanotate_minorest_version + ".")
	logger.info("phan_out version found is v" + str(phanotate_major_version) +"." + str(phanotate_minor_version) +"."+phanotate_minorest_version +".")

	if phanotate_major_version < 1:
		sys.exit("Phanotate is too old - please reinstall pharokka.")
	if phanotate_minor_version < 5:
		sys.exit("Phanotate is too old - please reinstall pharokka.")

	print("Phanotate version is ok.")
	logger.info("Phanotate version is ok.")

	#############
	# mmseqs
	#############
	try:
		process = sp.Popen(["mmseqs"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	except:
		sys.exit("MMseqs2 not found. Please reinstall pharokka.")
	mmseqs_out, _ = process.communicate()
	mmseqs_out = mmseqs_out.decode()

	version_line = []

	for line in mmseqs_out.split("\n"):
		if "Version" in line:
			version_line.append(line)

	mmseqs_version = version_line[0].split(' ')[2]
	mmseqs_major_version = int(mmseqs_version.split('.')[0])
	mmseqs_minor_version = int(mmseqs_version.split('.')[1])

	print("MMseqs2 version found is v" + str(mmseqs_major_version) +"." + str(mmseqs_minor_version) +".")
	logger.info("MMseqs2 version found is v" + str(mmseqs_major_version) +"." + str(mmseqs_minor_version) +".")

	if mmseqs_major_version != 13:
		sys.exit("MMseqs2 is the wrong version. Please install v13.45111")
	if mmseqs_minor_version != 45111:
		sys.exit("MMseqs2 is the wrong version. Please install v13.45111")

	print("MMseqs2 version is ok.")
	logger.info("MMseqs2 version is ok.")

	#############
	# trnascan
	#############
	try:
		process = sp.Popen(["tRNAscan-SE", "-h"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	except:
		sys.exit("tRNAscan-SE not found. Please reinstall pharokka.")

	trna_out, _ = process.communicate()
	trna_out = trna_out.decode()

	version_line = []

	for line in trna_out.split("\n"):
		if "tRNAscan-SE" in line:
			version_line.append(line)

	trna_version = version_line[0].split(' ')[1]
	trna_major_version = int(trna_version.split('.')[0])
	trna_minor_version = int(trna_version.split('.')[1])
	trna_minorest_version = int(trna_version.split('.')[2])

	print("tRNAscan-SE version found is v" + str(trna_major_version) +"." + str(trna_minor_version) +"." + str(trna_minorest_version) +".")
	logger.info("tRNAscan-SE version found is v" + str(trna_major_version) +"." + str(mmseqs_minor_version) +"."+ str(trna_minorest_version) +".")

	if trna_major_version != 2:
		sys.exit("tRNAscan-SE is the wrong version. Please re-install pharokka.")
	if trna_minor_version != 0:
		sys.exit("tRNAscan-SE is the wrong version. Please re-install pharokka.")
	if trna_minorest_version < 9:
		sys.exit("tRNAscan-SE is the wrong version. Please re-install pharokka.")

	print("tRNAscan-SE version is ok.")
	logger.info("tRNAscan-SE version is ok.")


	#############
	# minced
	#############
	try:
		process = sp.Popen(["minced", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	except:
		sys.exit("MinCED not found. Please reinstall pharokka.")

	minced_out, _ = process.communicate()
	minced_out = minced_out.decode()

	version_line = []

	for line in minced_out.split("\n"):
		if "minced" in line:
			version_line.append(line)

	minced_version = version_line[0].split(' ')[1]
	minced_major_version = int(minced_version.split('.')[0])
	minced_minor_version = int(minced_version.split('.')[1])
	minced_minorest_version = int(minced_version.split('.')[2])

	print("MinCED version found is v" + str(minced_major_version) +"." + str(minced_minor_version) +"." + str(minced_minorest_version) +".")
	logger.info("MinCED version found is v" + str(minced_major_version) +"." + str(minced_minor_version) +"."+ str(minced_minorest_version) +".")

	if minced_major_version != 0:
		sys.exit("MinCED is the wrong version. Please re-install pharokka.")
	if minced_minor_version != 4:
		sys.exit("MinCED is the wrong version. Please re-install pharokka.")
	if minced_minorest_version < 2:
		sys.exit("MinCED is the wrong version. Please re-install pharokka.")

	print("MinCED version is ok.")
	logger.info("MinCED version is ok.")


	#############
	# aragorn
	#############
	try:
		process = sp.Popen(["aragorn", "-h"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	except:
		sys.exit("ARAGORN not found. Please reinstall pharokka.")

	aragorn_out, _ = process.communicate()
	aragorn_out = aragorn_out.decode()

	version_line = []

	for line in aragorn_out.split("\n"):
		if "Dean Laslett" in line:
			version_line.append(line)

	aragorn_version = version_line[0].split(' ')[1]
	# strip off v
	aragorn_version = aragorn_version[1:]
	aragorn_major_version = int(aragorn_version.split('.')[0])
	aragorn_minor_version = int(aragorn_version.split('.')[1])
	aragorn_minorest_version = int(aragorn_version.split('.')[2])

	print("ARAGORN version found is v" + str(aragorn_major_version) +"." + str(aragorn_minor_version) +"." + str(aragorn_minorest_version) +".")
	logger.info("ARAGORN version found is v" + str(aragorn_major_version) +"." + str(aragorn_minor_version) +"."+ str(aragorn_minorest_version) +".")

	if aragorn_major_version != 1:
		sys.exit("ARAGORN is the wrong version. Please re-install pharokka.")
	if aragorn_minor_version != 2:
		sys.exit("ARAGORN is the wrong version. Please re-install pharokka.")
	if aragorn_minorest_version < 41:
		sys.exit("ARAGORN is the wrong version. Please re-install pharokka.")

	print("ARAGORN version is ok.")
	logger.info("ARAGORN version is ok.")



	#############
	# mash
	#############
	try:
		process = sp.Popen(["mash", "help"], stdout=sp.PIPE, stderr=sp.STDOUT) 
	except:
		sys.exit("mash not found. Please reinstall pharokka.")

	mash_out, _ = process.communicate()
	mash_out = mash_out.decode()

	version_line = []

	for line in mash_out.split("\n"):
		if "version" in line:
			version_line.append(line)

	mash_version = version_line[0].split(' ')[2]

	mash_major_version = int(mash_version.split('.')[0])
	mash_minor_version = int(mash_version.split('.')[1])

	print("mash version found is v" + str(mash_major_version) +"." + str(mash_minor_version) +"." )
	logger.info("mash version found is v" + str(mash_major_version) +"." + str(mash_minor_version) +".")

	if mash_major_version != 2:
		sys.exit("mash is the wrong version. Please re-install pharokka.")
	if mash_minor_version < 2:
		sys.exit("mash is the wrong version. Please re-install pharokka.")

	print("mash version is ok.")
	logger.info("mash version is ok.")