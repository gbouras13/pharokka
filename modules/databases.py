#!/usr/bin/env python3
import os
import sys
import subprocess as sp

def instantiate_dirs():
	out_dir = "output/"
	if os.path.isdir(out_dir) == False:
		os.mkdir(out_dir)
	db_dir = "databases/"
	if os.path.isdir(db_dir) == False:
		os.mkdir(db_dir)
	mmseqs_dir = "output/mmseqs/"
	if os.path.isdir(mmseqs_dir) == False:
		os.mkdir(mmseqs_dir)

def get_phrog_mmseqs():
    print("Getting PHROGs MMSeqs DB")
    db_dir = "databases/"
    filepath = "https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz"
    tarball = "phrogs_mmseqs_db.tar.gz"
    folder = "phrogs_mmseqs_db"
    
    # get tarball if not already present
    if os.path.isfile(os.path.join(db_dir,tarball)) == True: 
         print("PHROGs Database already downloaded")
        # download tarball and untar
    else:
        try:
            sp.call(["wget", filepath, "-P", db_dir])
        except:
            sys.stderr.write("Error: PHROGs Database not found - link likely broken\n")  
            return 0

    # delete folder if it exists already
    if os.path.isfile(os.path.join(db_dir,folder)) == True:
        sp.call(["rm", os.path.join(db_dir,folder)])
    
    # download untar -C for specifying the directory
    sp.call(["tar", "-xzf", os.path.join(db_dir, tarball), "-C", db_dir])


def get_phrog_annot_table():
    print("Getting PHROGs Annotation Table")
    db_dir = "databases/"
    filepath = "https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v3.tsv"
    file = "phrog_annot_v3.tsv"
    #if the file already exists
    if os.path.isfile(os.path.join(db_dir,file)) == True:
        print("PHROGs annotation file already downloaded")
    else:
        try:
            sp.call(["wget", filepath, "-P", db_dir])
        except:
            sys.stderr.write("Error: PHROGs annotation file not found - link likely broken\n")  
            return 0

