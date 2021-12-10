#!/usr/bin/env python3
"""The setup script."""

from modules import databases

databases.instantiate_dirs()
databases.get_phrog_mmseqs()
databases.get_phrog_annot_table()