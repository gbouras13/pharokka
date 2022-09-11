#!/usr/bin/env python3
import os
import sys
import subprocess as sp

def instantiate_install(db_dir):
    instantiate_dir(db_dir)
    get_phrog_mmseqs(db_dir)
    get_phrog_annot_table(db_dir)
    get_vfdb(db_dir)
    get_card(db_dir)


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
    filepath = "https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv"
    file = "phrog_annot_v4.tsv"
    #if the file already exists
    if os.path.isfile(os.path.join(db_dir,file)) == True:
        print("PHROGs annotation file already downloaded")
    else:
        try:
            sp.call(["curl", filepath, "-o", os.path.join(db_dir,file)])
        except:
            sys.stderr.write("Error: PHROGs annotation file not found - link likely broken\n")  
            return 0

def get_vfdb(db_dir):
    print("Getting VFDB Annotation Table")
    filepath = "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz"
    file = "VFDB_setB_pro.fas.gz"
    #if the file already exists
    if os.path.isfile(os.path.join(db_dir,"vfdb", "vfdb")) == True:
        print("VFDB already downloaded")
    else:
        try:
            instantiate_dir(os.path.join(db_dir, "vfdb"))
            sp.call(["curl", filepath, "-o", os.path.join(db_dir,"vfdb",file)])
            sp.Popen(["gunzip",  os.path.join(db_dir,"vfdb", file)], stdout=sp.PIPE)
            sp.call(["mmseqs", "createdb", os.path.join(db_dir, "vfdb", "VFDB_setB_pro.fas"), os.path.join(db_dir, "vfdb", "vfdb")])
        except:
            sys.stderr.write("Error: VFDB  not found - link likely broken\n")  
            return 0

def get_card(db_dir):
    print("Getting CARD Database Annotation Table")
    filepath = "https://card.mcmaster.ca/download/0/broadstreet-v3.2.4.tar.bz2"
    file = "card.tar.bz2"
    #if the file already exists
    if os.path.isfile( os.path.join(db_dir, "CARD_mmseqs", "CARD")) == True:
        print("CARD already downloaded")
    else:
        try:
            # make the CARD dir
            instantiate_dir(os.path.join(db_dir, "CARD"))
            instantiate_dir(os.path.join(db_dir, "CARD_mmseqs"))
            # download the database 
            sp.call(["curl", filepath, "-o", os.path.join(db_dir,"CARD",file)])
            # untar 
            sp.call(["tar", "-xzf", os.path.join(db_dir,"CARD",file), "-C",os.path.join(db_dir,"CARD") ])
            # create mmseqs db
            sp.call(["mmseqs", "createdb", os.path.join(db_dir, "CARD", "protein_fasta_protein_homolog_model.fasta"), os.path.join(db_dir, "CARD_mmseqs", "CARD")])
        except:
            sys.stderr.write("Error: CARD  not found - link likely broken\n")  
            return 0