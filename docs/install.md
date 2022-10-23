# Installation

**pharokka v1.1.0 is now available on bioconda**

**Important: MMseqs2 has recently been updated to v14-7e284 and that version will not work with Pharokka. Please read below**

* The MMseqs2 team have recently changed the internal MMseqs2 profile format in the new MMseqs2 version v14-7e284.
* This means that the pharokka database will not work with MMseqs2 v14-7e284.
* As a result, pharokka needs to be run with MMseqs2 v13.4511.
* As of v1.1.0, the pharokka bioconda dependencies are fixed to ensure that MMseqs2 is v13.4511. 
* However, if you run into issues with pharokka (particularly pre v1.1.0 versions), please check that MMseqs2 is v13.4511 is installed (it will be clear in the pharokka***.log output file). I would recommend a fresh pharokka install in this instance, or trying `conda install -c bioconda pharokka mmseqs2==13.4511`.
* If you are installing pharokka from the git repository, the environment.yml file has been changed to fix MMseqs2 v13.4511, so proceed as below. 

The easiest way to install pharokka is via conda. For inexperienced command line users, this method is highly recommended.

`conda install -c bioconda pharokka`

This will install all the dependencies along with pharokka. The dependencies are listed in environment.yml.

If conda is taking a long time to solve the environment, try using mamba:

```
conda install mamba
mamba install -c bioconda pharokka
```

Alternatively, the development version of pharokka can be installed manually via github. 

`git clone https://github.com/gbouras13/pharokka.git`

The dependencies found in environment.yml will then need to be installed manually.

For example using conda to install the required dependencies:

```
git clone https://github.com/gbouras13/pharokka.git
cd pharokka
conda env create -f environment.yml
conda activate pharokka_env
```

And then to run pharokka (assuming you are still in the pharokka directory)

```
./bin/install_databases.py -h
./bin/pharokka.py -h
```

# Database Installation

* v1.0.0 onwards adds VFDB (current as of 15-09-22) and CARD (v3.2.4) databases for virulence factor and AMR gene identification.
* These should install using the install_databases.py script, with the databases downloaded from a Zenodo repository.
* You will need to re-install the databases if you updating from an earlier version of pharokka than v1.0.0. The database should work for all versions from v1.0.0 and afterwards.
* If the script does not work, you an alternatively download the databases manually from Zenodo at https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz and untar the directory in a location of your choice. Please see the Installation Section for more details.

To install the pharokka database to the default directory:

`install_databases.py -d`

If you would like to specify a different database directory (recommended), that can be achieved as follows:

`install_databases.py -o <path/to/databse_dir>`

If this does not work, you an alternatively download the databases from Zenodo at https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz and untar the directory in a location of your choice.

If you prefer to use the command line:

```
wget "https://zenodo.org/record/7081772/files/pharokka_database_v1.0.0.tar.gz"
tar -xzf pharokka_database_v1.0.0.tar.gz
```

which will create a directory called "pharokka_database_v1.0.0" containing the databases.

# Beginner Conda Installation

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

5. After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

`conda install mamba`

 6. Finally, I would recommend installing pharokka into a fresh environment. For example to create an environment called pharokkaENV with pharokka installed:

```
mamba create -n pharokkaENV pharokka
conda activate pharokkaENV
install_databases.py -h
pharokka.py -h
```
