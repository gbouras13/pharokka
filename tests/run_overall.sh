#!/bin/bash

# install_databases.py -o tests/test_data/database

test_data="tests/test_data/overall"
db_dir="tests/test_data/database"
custom_dir="tests/test_data/custom_db"
out_dir="output"
mkdir -p $out_dir

# normal
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1 -t 8 -f

# AMR
pharokka.py -i $test_data/AMR_example/NC_007458.fasta -d $db_dir -o $out_dir/NC_007458 -t 8 -f

# genbank
pharokka.py -i $test_data/Standard_examples/SAOMS1_Output/pharokka.gbk -d $db_dir -o $out_dir/SAOMS1_genbank -t 8 -f --genbank

# meta genbank
pharokka.py -i $test_data/genbank_examples/hundred_microviruses.gbk -d $db_dir -o $out_dir/hundred_microviruses_genbank -m -t 8 -f --genbank

# CRISPR
pharokka.py -i $test_data/CRISPR_example/Biggiephage_A_fullcontig_CasΦ1.fasta -d $db_dir -o $out_dir/Biggiephage -t 8 -f

# tmRNA
pharokka.py -i $test_data/tmRNA_example/NC_051700.fasta -d $db_dir -o $out_dir/tmRNA -t 8 -f

# vfdb
pharokka.py -i $test_data/VFDB_example/NC_004617.fasta -d $db_dir -o $out_dir/VFDB -t 8 -f

# normal prod
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_prod -g prodigal -t 8 -f

# stop recode
pharokka.py -i $test_data/stop_recoding/table_4/SRR1747055_scaffold_7.fa -d $db_dir -o $out_dir/table4 -g prodigal -t 8 -f -c 4

# normal
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_dnaap -t 8 -f --dnaapler

# fast
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_fast -t 8 -f --fast

# mmseqs2_only
pharokka.py -i $test_data/Standard_examples/SAOMS1.fasta -d $db_dir -o $out_dir/SAOMS1_mmseqs2_only -t 8 -f --mmseqs2_only

# custom db
pharokka.py -i $test_data/custom_examples/MH649026.fasta -d $db_dir -o $out_dir/MH649026 -t 8 -f --custom_hmm $custom_dir/microvirus.h3m

# custom db meta
pharokka.py -i $test_data/custom_examples/hundred_microviruses.fasta -d $db_dir -o $out_dir/hundred_microviruses -t 8 -m -s -f --custom_hmm $custom_dir/microvirus.h3m --fast

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

pharokka_proteins.py -i tests/test_data/proteins/phanotate.faa -d $db_dir -o $out_dir/proteins_test -t 8 -f
pharokka_proteins.py -i tests/test_data/proteins/phanotate.faa -d $db_dir -o $out_dir/proteins_test_mmseqs2_only -t 8 -f --mmseqs2_only

