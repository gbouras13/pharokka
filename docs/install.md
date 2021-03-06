**pharokka v0.1.7 is now available on bioconda**

The easiest way to install pharokka is via conda.

`conda install -c bioconda pharokka`

This will install all the dependencies along with pharokka.

Alternatively, pharokka can be installed manually via github

`git clone https://github.com/gbouras13/pharokka.git`

The dependencies found in requirements.yml will then need to be installed manually.

For example using conda:

```
git clone https://github.com/gbouras13/pharokka.git
cd pharokka
conda env create -f environment.yml
conda activate pharokka_env
install_databses.py -h
pharokka.py -h
```

Before running pharokka, the PHROGs database needs to be downloaded using

`install_databases.py -d `

If you would like to specify a different database directory (recommended), that can be achieved as follows:

`install_databases.py -o <path/to/databse_dir>`

If you have trouble downloading the databases using install_databases.py, they can be manually downloaded from the PHROGs website and placed in a directory of your choice:

* [https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz](https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz)
* [https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv](https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv).
