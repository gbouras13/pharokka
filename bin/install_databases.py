#!/usr/bin/env python3
"""The setup script."""

from modules import databases
import os


db_dir = os.path.join(os.path.dirname(__file__),'../',"databases/")  

databases.instantiate_install(db_dir)