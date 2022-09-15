**pharokka v1.0.0 is now available on bioconda**

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


Beginner Conda Installation
--------

If you are new to using the command-line, please install conda using the following instructions.

1. Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

3. Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`

4. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

5. After this, conda should be installed (you may need to restart your terminal). I would recommend installing pharokka into a fresh environment e.g. to create an environment called pharokkaENV with pharokka installed:

```
conda create -n pharokkaENV pharokka
conda activate pharokkaENV
```

Database Installation
-----------------

Before running pharokka, the PHROGs database needs to be downloaded using

`install_databases.py -d`

If you would like to specify a different database directory (recommended), that can be achieved as follows:

`install_databases.py -o <path/to/databse_dir>`

v1.0.0 adds VFDB and CARD databases for virulence factor and AMR gene identification. These should install using the install_databases.py script as outlined above. You will need to run this before running pharokka v1.0.0.

If this does not work, you an alternatively download the databases from Zenodo at https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz and unzip the directory in a location of your choice.

If you prefer to use the command line:

```
wget "https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz"
tar -xzf pharokka_database_v1.0.0.tar.gz
```

which will create a directory called "pharokka_database_v1.0.0" containing the databases.

If you have already downloaded databases for earlier versions of pharokka, these will need to be re-downloaded.
