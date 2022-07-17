
pharokka
===============

Fast Phage Annotation Program
------------

pharokka is designed for rapid standardised annotation of bacteriophages.

If you are looking for rapid standardised annotation of prokaryotes, please use prokka (https://github.com/tseemann/prokka), which inspired the creation of pharokka.

Method
----
Briefly, default gene prediction is done using PHANOTATE (https://github.com/deprekate/PHANOTATE) and function annotation is based on the PHROGs database (https://phrogs.lmge.uca.fr).

The main output is a gff file that is suitable for use downstream pangenomic pipelines such as Roary (https://sanger-pathogens.github.io/Roary/).

The other important output is `cds_functions.tsv`, which includes counts of CDSs, tRNAs, and functions assigned to CDSs according to the PHROGs database.

For full documentation, please visit https://pharokka.readthedocs.io.

Usage
------
The easiest way to install pharokka is via conda using

`conda install pharokka -c gbouras13`

This will install all the dependencies along with pharokka.

Alternatively, the development version of pharokka can be installed manually via github.

`git clone https://github.com/gbouras13/pharokka.git`

The dependencies found in environment.yml will then need to be installed manually.

For example using conda:

```
conda env create -f environment.yml
conda activate pharokka_env
git clone https://github.com/gbouras13/pharokka.git
cd pharokka
```

Running pharokka
--------

First the PHROGs databases need to be installed

`install_databses.py -d Y`

If you would like to specify a different database directory (recommended), that can be achieved as follows:

`install_databases.py -d N -o "<path/to/databse_dir>`

If you have trouble downloading the databases using `install_databases.py`, they can be manually downloaded from the PHROGs website links, untared and placed in a directory of your choice:
* https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz
* https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz
* https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv.

Once the databases have finished downloading, run pharokka

`pharokka.py -i <fasta file> -o <output folder> -t <threads>`

To specify a prefix for the output files:

`pharokka.py -i <fasta file> -o <output folder> -t <threads> -p <prefix>`

To specify a different database directory (recommended):

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads> -p <prefix>`

To overwrite an existing output directory, use -f

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -p <prefix> -f`

To use Prodigal instead of PHANOTATE use `-g prodigal`

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal`

pharokka should work with metagenome assembled viral contigs with PHANOTATE automatically. With prodigal, please add the -m flag

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/database_dir> -t <threads>  -g prodigal -m`

pharokka defaults to 1 thread.

System
------
pharokka has been tested on Linux and MacOS (M1 and Intel).

Time
--------
On a standard 16GB laptop specifying 8 threads, pharokka should take between 5-20 minutes, depending on the genome size.

Bugs and Suggestions
--------
If you come across bugs with pharokka, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au

Citation
--------
If you use pharokka, please also cite:

* PHANOTATE (McNair K., Zhou C., Dinsdale E.A., Souza B., Edwards R.A. (2019) "PHANOTATE: a novel approach to gene identification in phage genomes", Bioinformatics, https://doi.org/10.1093/bioinformatics/btz26).
* tRNAscan-SE: Chan, P.P., Lin, B.Y., Mak, A.J. and Lowe, T.M. (2021) "tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes", Nucleic Acids Res., https://doi.org/10.1093/nar/gkab688.
* Steinegger M. and Soeding J. (2017), "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets", Nature Biotechnology (https://doi.org/10.1038/nbt.3988).
* Terzian P., Olo Ndela E., Galiez C., Lossouarn J., Pérez Bucio R.E., Mom R., Toussaint A., Petit M.A., Enault F., "PHROG : families of prokaryotic virus proteins clustered using remote homology", NAR Genomics and Bioinformatics, (2021), (https://doi.org/10.1093/nargab/lqab067).
