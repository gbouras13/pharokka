#!/usr/bin/env python3
import os
import sys
import subprocess as sp

PHROG_DB_NAMES = ['phrogs_db','phrogs_db.dbtype',
'phrogs_db.index',
'phrogs_profile_db',
'phrogs_profile_db.dbtype',
'phrogs_profile_db.index',
'phrogs_profile_db_consensus',
'phrogs_profile_db_consensus.dbtype',
'phrogs_profile_db_consensus.index',
'phrogs_profile_db_h',
'phrogs_profile_db_h.index',
'phrogs_profile_db_seq',
'phrogs_profile_db_seq.dbtype',
'phrogs_profile_db_seq.index',
'phrogs_profile_db_seq_h',
'phrogs_profile_db_seq_h.index']

VFDB_DB_NAMES = ['VFDB_setB_pro.fas'
'vfdb',
'vfdb.dbtype',
'vfdb.index',
'vfdb.lookup',
'vfdb.source',
'vfdb_h',
'vfdb_h.dbtype',
'vfdb_h.index']

CARD_DB_NAMES = [
'CARD',
'CARD.dbtype',
'CARD.index',
'CARD.lookup',
'CARD.source',
'CARD_h',
'CARD_h.dbtype',
'CARD_h.index']

def instantiate_install(db_dir):
    instantiate_dir(db_dir)
    downloaded_flag = check_db_installation(db_dir)
    if downloaded_flag == True:
        print("All Databases have been Downloaded and Checked")
    else:
        get_database_zenodo(db_dir)

def instantiate_dir(db_dir):
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)
 
def check_db_installation(db_dir):

    downloaded_flag = True

    for file_name in PHROG_DB_NAMES:
        path = os.path.isfile(os.path.join(db_dir,'phrogs_mmseqs_db', file_name))
        if os.path.isfile(path) == False:
            print("PHROGs Databases Need to be Downloaded")
            downloaded_flag = False
            break
    # VFDB
    for file_name in VFDB_DB_NAMES:
        path = os.path.isfile(os.path.join(db_dir,'vfdb', file_name))
        if os.path.isfile(path) == False:
            print("VFDB Databases Need to be Downloaded")
            downloaded_flag = False
            break
    # CARD
    for file_name in CARD_DB_NAMES:
        path = os.path.isfile(os.path.join(db_dir,'CARD_mmseqs', file_name))
        if os.path.isfile(path) == False:
            print("CARD Databases Need to be Downloaded")
            downloaded_flag = False
            break
    # annot.tsv
    path = os.path.isfile(os.path.join(db_dir,'phrog_annot_v4.tsv'))
    if os.path.isfile(path) == False:
            print("PHROGs Annotation File Needs to be Downloaded")
            downloaded_flag = False
    
    return downloaded_flag
    

def get_database_zenodo(db_dir):
    print("Downloading Pharokka Database")
    file = 'pharokka_v0.1.11_databases.zip'
    url = "https://zenodo.org/record/7080544/files/pharokka_v0.1.11_databases.zip"
    try:
        sp.call(["rm", os.path.join(db_dir,file)])
        sp.call(["curl", url, "-o", os.path.join(db_dir,file)])
        sp.call(["unzip", '-j', os.path.join(db_dir,file), '-d', db_dir])
    except:
        sys.stderr.write("Error: Pharokka Database Install Failed. \n Please try again or use the manual option detailed at https://github.com/gbouras13/pharokka.git \n")  
        return 0


