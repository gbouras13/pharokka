The easiest way to install phrokka is via conda using

`conda install phrokka -c gbouras13`

This will install all the dependencies along with phrokka.

Alternatively, phrokka can be installed manually via github

`git clone https://github.com/gbouras13/phrokka.git`

The dependencies found in other_files/requirements.yaml will then need to be installed manually.

For example using conda:

```
conda env create -f environment.yml -n phrokka_env
conda activate phrokka_env
git clone https://github.com/gbouras13/phrokka.git
```

Before running phrokka, the PHROGs database needs to be downloaded using

`install_databases.py -d Y`


If you would like to specify a different database directory, that can be achieved as follows:

`install_databases.py -d N -o "<path/to/databse_dir>`
