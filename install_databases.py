#!/usr/bin/env python3
"""The setup script."""

from modules import databases
import global_variables

db_dir = global_variables.DATABASE_DIR

databases.instantiate_install(db_dir)