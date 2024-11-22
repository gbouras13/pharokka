# Installation

## Conda

The easiest way to install `pharokka` is via conda. For inexperienced command line users, this method is highly recommended.

```bash
conda install -c bioconda pharokka
```

This will install all the dependencies along with `pharokka`. The dependencies are listed in environment.yml.

If you have a MacOS system with M1/M2/M3 Apple Silicon, try this 

```bash
conda create --platform osx-64 --name pharokkaENV -c bioconda pharokka
conda activate pharokkaENV
```

## Pip

As of v1.4.0, you can also install the python components of `pharokka` with pip.

```bash
pip install pharokka
```

You will still need to install the non-python dependencies manually.

## Container

If you have Docker/Singularity/Apptainer installed, you can use the [biocontainers container](https://quay.io/repository/biocontainers/pharokka?tab=tags) (yes, every bioconda package has one!)

You might find this useful if you have trouble with conda environments.

For example to install `pharokka v1.7.3` with Singularity:

```
IMAGE_DIR="<the directory you want the .sif file to be in >"
# e.g. to pull into the working directory
IMAGE_DIR=$PWD
singularity pull --dir $IMAGE_DIR docker://quay.io/biocontainers/pharokka:1.7.3--pyhdfd78af_0
```

* Then to run, you use the same commands but prepended with `singularity exec <.sif file>` e.g.:

```
containerImage="$IMAGE_DIR/pharokka_1.7.3--pyhdfd78af_0.sif"
singularity exec $containerImage pharokka.py -h
```

## Source

Alternatively, the development version of `pharokka` (which may include new, untested features) can be installed manually via github. 

```bash
git clone https://github.com/gbouras13/pharokka.git
cd pharokka
pip install -e .
pharokka.py --help
```

The dependencies found in environment.yml will then need to be installed manually.

For example using conda to install the required dependencies:

```bash
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

If this does not work, you an alternatively download the databases from Zenodo at https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz and untar the directory in a location of your choice.

If you prefer to use the command line:

```bash
wget "https://zenodo.org/record/8276347/files/pharokka_v1.4.0_databases.tar.gz"
tar -xzf pharokka_v1.4.0_databases.tar.gz
```

which will create a directory called "pharokka_v1.4.0_databases" containing the databases.

# Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

1. Install Conda - I would recommend [miniforge](https://github.com/conda-forge/miniforge).
2. Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniforge](https://github.com/conda-forge/miniforge) repository).

`wget https://github.com/conda-forge/miniforge/releases/download/24.9.2-0/Miniforge3-24.9.2-0-Linux-x86_64.sh`

For Mac Intel:

`wget https://github.com/conda-forge/miniforge/releases/download/24.9.2-0/Miniforge3-24.9.2-0-MacOSX-x86_64.sh`

For Mac M1/M2/M3/M4

`wget https://github.com/conda-forge/miniforge/releases/download/24.9.2-0/Miniforge3-24.9.2-0-MacOSX-arm64.sh`

3. Install miniforge and follow the prompts.

`sh Miniforge3-24.9.2-0-Linux-x86_64.sh`

4. After installation is complete, you should add the following channels to your conda configuration:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

 5. Finally, I would recommend installing pharokka into a fresh environment. For example to create an environment called pharokkaENV with pharokka installed:

```bash
conda create -n pharokkaENV pharokka
conda activate pharokkaENV
install_databases.py -h
pharokka.py -h
```

If you have a Mac with Apple Silicon (M1-M4), try

```bash
conda create --platform osx-64 -n pharokkaENV pharokka
conda activate pharokkaENV
install_databases.py -h
pharokka.py -h
```