#!/usr/bin/env python3
import os
import sys
import subprocess as sp

def instantiate_install(db_dir):
    instantiate_dir(db_dir)
    get_phrog_mmseqs(db_dir)
    get_phrog_annot_table(db_dir)
    get_phrog_hhmer(db_dir)

def instantiate_dir(db_dir):
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)
 
def get_phrog_mmseqs(db_dir):
    print("Getting PHROGs MMSeqs DB")
    filepath = "https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz"
    tarball = "phrogs_mmseqs_db.tar.gz"
    folder = "phrogs_mmseqs_db"
    
    # get tarball if not already present
    if os.path.isfile(os.path.join(db_dir,tarball)) == True: 
         print("PHROGs Database already downloaded")
        # download tarball and untar
    else:
        try:
            sp.call(["curl", filepath, "-o", os.path.join(db_dir,tarball)])
        except:
            sys.stderr.write("Error: PHROGs MMSeqs Database not found - link likely broken\n")  
            return 0

    # delete folder if it exists already
    if os.path.isfile(os.path.join(db_dir,folder)) == True:
        sp.call(["rm", os.path.join(db_dir,folder)])
    
    # download untar -C for specifying the directory
    sp.call(["tar", "-xzf", os.path.join(db_dir, tarball), "-C", db_dir])


def get_phrog_annot_table(db_dir):
    print("Getting PHROGs Annotation Table")
    filepath = "https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v3.tsv"
    file = "phrog_annot_v3.tsv"
    #if the file already exists
    if os.path.isfile(os.path.join(db_dir,file)) == True:
        print("PHROGs annotation file already downloaded")
    else:
        try:
            sp.call(["curl", filepath, "-o", os.path.join(db_dir,file)])
        except:
            sys.stderr.write("Error: PHROGs annotation file not found - link likely broken\n")  
            return 0

def get_phrog_hhmer(db_dir):
    print("Getting PHROGs HHmer DB")
    filepath = "https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz"
    tarball = "phrogs_hhsuite_db.tar.gz"
    folder = "phrogs_hhsuite_db"
    
    # get tarball if not already present
    if os.path.isfile(os.path.join(db_dir,tarball)) == True: 
         print("PHROGs Database already downloaded")
        # download tarball and untar
    else:
        try:
            sp.call(["curl", filepath, "-o", os.path.join(db_dir,tarball)])
        except:
            sys.stderr.write("Error: PHROGs HMMer Database not found - link likely broken\n")  
            return 0

    # delete folder if it exists already
    if os.path.isfile(os.path.join(db_dir,folder)) == True:
        sp.call(["rm", os.path.join(db_dir,folder)])
    
    # download untar -C for specifying the directory
    sp.call(["tar", "-xzf", os.path.join(db_dir, tarball), "-C", db_dir])
