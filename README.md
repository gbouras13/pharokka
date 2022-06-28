
pharokka
===============

Fast Phage Annotation Pipeline
------------

pharokka is designed for rapid, standard annotation of bacteriophages based on the PHROGs database (https://phrogs.lmge.uca.fr).

For full documentation, please visit https://pharokka.readthedocs.io.

The easiest way to install pharokka is via conda using

`conda install pharokka -c gbouras13`

This will install all the dependencies along with pharokka.

Alternatively, pharokka can be installed manually via github.

`git clone https://github.com/gbouras13/pharokka.git`

The dependencies found in environment.yml will then need to be installed manually.

For example using conda:

```
cd pharokka
conda env create -f environment.yml
conda activate pharokka_env
git clone https://github.com/gbouras13/pharokka.git
```

Running pharokka
--------

First the PHROGs databases need to be installed

`install_databses.py -d Y`

If you would like to specify a different database directory, that can be achieved as follows:

`install_databases.py -d N -o "<path/to/databse_dir>`

If you have trouble downloading the databases using install_databases.py, they can be manually downloaded from the PHROGs website and placed in a directory of your choice https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v3.tsv.

Once the databases have finished downloading, run pharokka

`pharokka.py -i <fasta file> -o <output folder> -t <threads>`

To specify a different database directory:

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads> `

To overwrite an existing output directory, use -f

`pharokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir> -t <threads>  -f`

pharokka defaults to 1 thread.
