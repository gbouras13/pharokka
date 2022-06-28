The easiest way to install pharokka is via conda using

`conda install pharokka -c gbouras13`

This will install all the dependencies along with pharokka.

Alternatively, pharokka can be installed manually via github

`git clone https://github.com/gbouras13/pharokka.git`

The dependencies found in other_files/requirements.yaml will then need to be installed manually.

For example using conda:

```
cd pharokka
conda env create -f environment.yml
conda activate pharokka_env
git clone https://github.com/gbouras13/pharokka.git
```

Before running pharokka, the PHROGs database needs to be downloaded using

`install_databases.py -d Y`


If you would like to specify a different database directory, that can be achieved as follows:

`install_databases.py -d N -o "<path/to/databse_dir>`

If you have trouble downloading the databases using install_databases.py, they can be manually downloaded from the PHROGs website and placed in a directory of your choice https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v3.tsv.
