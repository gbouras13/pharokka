# Installation

## Conda

The easiest way to install `pharokka` is via conda. For inexperienced command line users, this method is highly recommended.

```
conda install -c bioconda pharokka
```

This will install all the dependencies along with `pharokka`. The dependencies are listed in environment.yml.

If conda is taking a long time to solve the environment, try using mamba:

```
conda install mamba
mamba install -c bioconda pharokka
```

## Pip

As of v1.4.0, you can also install the python components of `pharokka` with pip.

```
pip install pharokka
```

You will still need to install the non-python dependencies manually.

## Source

Alternatively, the development version of `pharokka` (which may include new, untested features) can be installed manually via github. 

```
git clone https://github.com/gbouras13/pharokka.git
cd pharokka
pip install -e .
pharokka.py --help
```

The dependencies found in environment.yml will then need to be installed manually.

For example using conda to install the required dependencies:

```
conda env create -f environment.yml
conda activate pharokka_env
# assuming you are in the pharokka directory 
# installs pharokka from source
pip install -e .
pharokka.py --help
```

# Database Installation

* **Note v 1.4.0 implements a new database with PHROGs HMM profiles. You will need to update the Pharokka database to use v1.4.0 and higher**

To install the pharokka database to the default directory:

`install_databases.py -d`

If you would like to specify a different database directory (recommended), that can be achieved as follows:

`install_databases.py -o <path/to/databse_dir>`

If this does not work, you an alternatively download the databases from Zenodo at https://zenodo.org/record/8267900/files/pharokka_v1.4.0_databases.tar.gz and untar the directory in a location of your choice.

If you prefer to use the command line:

```
wget "https://zenodo.org/record/8267900/files/pharokka_v1.4.0_databases.tar.gz"
tar -xzf pharokka_v1.4.0_databases.tar.gz
```

which will create a directory called "pharokka_v1.4.0_databases" containing the databases.

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
