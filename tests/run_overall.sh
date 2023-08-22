#!/bin/bash

# install_databases.py -o tests/test_data/database

test_data="tests/test_data/overall"
db_dir="tests/test_data/database"
out_dir="output"
mkdir -p $out_dir

# normal
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir-o $out_dir/SAOMS1 -t 8 -f

# AMR
pharokka.py -i $test_data/AMR_example/NC_007458.fasta -d $db_dir -o $out_dir/NC_007458 -t 8 -f

# CRISPR
pharokka.py -i $test_data/CRISPR_example/Biggiephage_A_fullcontig_CasΦ1.fasta -d $db_dir -o $out_dir/Biggiephage -t 8 -f

# tmRNA
pharokka.py -i $test_data/tmRNA_example/NC_051700.fasta -d $db_dir -o $out_dir/tmRNA -t 8 -f

# tmRNA
pharokka.py -i $test_data/VFDB_example/NC_004617.fasta -d $db_dir -o $out_dir/VFDB -t 8 -f

# normal prod
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_prod -g prodigal -t 8 -f

# stop recode
pharokka.py -i $test_data/table_4/SRR1747055_scaffold_7.fa -d $db_dir -o $out_dir/SAOMS1_prod_table4 -g prodigal -t 8 -f -c 4

# normal
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_dnaap -t 8 -f --dnaapler

# fast
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_fast -t 8 -f --fast

# mmseqs2_only
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_mmseqs2_only -t 8 -f --mmseqs2_only

# meta
pharokka.py -i $test_data/Meta_example/fake_meta.fa -d $db_dir -o $out_dir/fake_meta -t 8 -f -m

# meta split
pharokka.py -i $test_data/Meta_example/fake_meta.fa -d $db_dir -o $out_dir/SAOMS1_split -t 8 -f -m -s

# bad meta
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_bad_meta -t 8 -f -m 

# terminase
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_terminase -t 8 -f --terminase --terminase_start 340 --terminase_strand neg

# mmseqs2_only
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_mmseqs2_only -t 8 -f --mmseqs2_only


# proteins

pharokka_proteins.py -i $test_data/proteins/phanotate.faa -d $db_dir -o $out_dir/proteins_test_mmseqs2_only -t 1 -f --mmseqs2_only