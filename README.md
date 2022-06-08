
phrokka
===============


Fast Phage Annotation Pipeline

For full documentation, please visit https://phrokka.readthedocs.io.

The easiest way to install phrokka is via conda using

`conda install phrokka -c gbouras13`

This will install all the dependencies along with phrokka.

Alternatively, phrokka can be installed manually via github.

`git clone https://github.com/gbouras13/phrokka.git`

The dependencies found in environment.yaml will then need to be installed manually.

For example using conda:

```
conda env create -f environment.yml
conda activate phrokka_env
git clone https://github.com/gbouras13/phrokka.git
```

Running phrokka
--------

First the PHROGs databases need to be installed

`install_databses.py -d Y`

If you would like to specify a different database directory, that can be achieved as follows:

`install_databases.py -d N -o "<path/to/databse_dir>`

Once the databases have finished downloading, run phrokka

`phrokka.py -i <fasta file> -o <output folder>`

To specify a different database directory:

`phrokka.py -i <fasta file> -o <output folder> -d <path/to/databse_dir>`

To overwrite an existing output directory, use -f

`phrokka.py -i <fasta file> -o <output folder> -f`
