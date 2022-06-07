import argparse
import os
from argparse import RawTextHelpFormatter



from modules import databases
import os

def get_db_input():
	usage = 'install_databases.py -d N -o "path/to/dir'
	parser = argparse.ArgumentParser(description='script to download required phrog databases', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d', '--default', action="store", help='Must be "Y" or "N". Determines whether you want databases stored in the default location, or in a custom directory',
    required=True)
	parser.add_argument('-o', '--outdir', action="store", help='Database Directory - will be created and must be specificed in conjunction with -d N')
	args = parser.parse_args()
	return args


args = get_db_input()
print(args.outdir)

if args.default == 'Y':
    db_dir = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    databases.instantiate_install(db_dir)
elif args.default == 'N':
    if args.outdir == None:
        print("-o was not specified, but -d N was.")
        print("Downloading databases to the default directory anyway")
        db_dir = os.path.join(os.path.dirname(__file__),'../',"databases/")  
    else:
        db_dir = args.outdir
    databases.instantiate_install(db_dir)
else:
    print("-d was not specified as Y or N. Please try again.")
    exit()
