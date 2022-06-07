#!/usr/bin/env python3
"""The setup script."""

from modules import databases
import os
import argparse
from argparse import RawTextHelpFormatter

def get_db_input():
	usage = 'install_databases.py -d N -o "path/to/dir'
	parser = argparse.ArgumentParser(description='script to download required phrog databases', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', '--default', action="store", help='Must be "Y" or "N". Determines whether you want databases stored in the default location, or in a custom directory',required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Database Directory - will be created and must be specificed in conjunction with -d N' )
	args = parser.parse_args()
	return args


args = get_db_input()

if args.default == 'Y':
    db_dir = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    databases.instantiate_install(db_dir)
elif args.default == 'N':
    db_dir = args.outdir
    databases.instantiate_install(db_dir)
else:
    print("-d was not specified as Y or N. Please try again.")
    exit()
