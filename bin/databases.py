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

VFDB_DB_NAMES = ['VFDB_setB_pro.fas',
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
        print("All Databases have already been Downloaded and Checked")
    else:
        get_database_zenodo(db_dir)

def instantiate_dir(db_dir):
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)
 
def check_db_installation(db_dir):

    downloaded_flag = True
    # PHROGS files
    for file_name in PHROG_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            print("PHROGs Databases Need to be Downloaded")
            downloaded_flag = False
            break
    # VFDB
    for file_name in VFDB_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            print("VFDB Databases Need to be Downloaded")
            downloaded_flag = False
            break
    # CARD
    for file_name in CARD_DB_NAMES:
        path = os.path.join(db_dir, file_name)
        if os.path.isfile(path) == False:
            print("CARD Databases Need to be Downloaded")
            downloaded_flag = False
            break
    # annot.tsv
    path = os.path.join(db_dir,'phrog_annot_v4.tsv')
    if os.path.isfile(path) == False:
            print("PHROGs Annotation File Needs to be Downloaded")
            downloaded_flag = False
    
    return downloaded_flag
    

def get_database_zenodo(db_dir):
    print("Downloading Pharokka Database")
    tarball = 'pharokka_v_1.0.0_databases.tar.gz'
    url = "https://zenodo.org/record/7080911/files/pharokka_v_1.0.0_databases.tar.gz"
    try:
        # remvoe the directory
        sp.call(["rm", "-rf", os.path.join(db_dir)])
        # make db dir
        sp.call(["mkdir", "-p", os.path.join(db_dir)])
        # download the tarball
        sp.call(["curl", url, "-o", os.path.join(db_dir,tarball)])
        # untar tarball into database directory
        sp.call(["tar", "-xzf", os.path.join(db_dir, tarball), "-C", db_dir, "--strip-components=1"])
        # remove tarball
        sp.call(["rm","-f", os.path.join(db_dir,tarball)])
    except:
        sys.stderr.write("Error: Pharokka Database Install Failed. \n Please try again or use the manual option detailed at https://github.com/gbouras13/pharokka.git \n")  
        return 0
