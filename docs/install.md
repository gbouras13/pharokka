# Installation

## Conda

The easiest way to install `pharokka` is via conda. For inexperienced command line users, this method is highly recommended.

```bash
conda install -c bioconda pharokka
```

This will install all the dependencies along with `pharokka`. The dependencies are listed in environment.yml.


## Pip

As of v1.4.0, you can also install the python components of `pharokka` with pip.

```bash
pip install pharokka
```

You will still need to install the non-python dependencies manually.

## Container

If you have Docker/Singularity/Apptainer installed, you can use the [biocontainers container](https://quay.io/repository/biocontainers/pharokka?tab=tags) (yes, every bioconda package has one!)

You might find this useful if you have trouble with conda environments.

For example to install `pharokka v1.9.1` with Singularity:

```
IMAGE_DIR="<the directory you want the .sif file to be in >"
# e.g. to pull into the working directory
IMAGE_DIR=$PWD
singularity pull --dir $IMAGE_DIR docker://quay.io/biocontainers/pharokka:1.9.1--pyhdfd78af_1
```

* Then to run, you use the same commands but prepended with `singularity exec <.sif file>` e.g.:

```
containerImage="$IMAGE_DIR/pharokka_1.9.1--pyhdfd78af_1.sif"
singularity exec $containerImage pharokka -h
```

## Source

Alternatively, the development version of `pharokka` (which may include new, untested features) can be installed manually via github. 

```bash
git clone https://github.com/gbouras13/pharokka.git
cd pharokka
pip install -e .
pharokka --help
```

The dependencies found in environment.yml will then need to be installed manually.

For example using conda to install the required dependencies:

```bash
conda env create -f environment.yml
conda activate pharokka_env
# assuming you are in the pharokka directory 
# installs pharokka from source
pip install -e .
pharokka --help
```

# Database Installation

* **Note: v1.8.0 uses a new MMseqs2 PHROG profile database format that is incompatible with the v1.4.0 database. If upgrading from v1.7.x or earlier, you must re-run `pharokka install` to fetch the updated database.**

To install the pharokka database to the default directory:

```bash
pharokka install -d
```

If you would like to specify a different database directory (recommended), that can be achieved as follows:

```bash
pharokka install -o <path/to/database_dir>
```

If this does not work, you can alternatively download the databases manually from Zenodo:

```bash
wget "https://zenodo.org/record/17110353/files/pharokka_v1.8.0_databases.tar.gz"
tar -xzf pharokka_v1.8.0_databases.tar.gz
```

This will create a directory called `pharokka_v1.8.0_databases` containing the databases.

# Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

1. Install Conda - I would recommend [miniforge](https://github.com/conda-forge/miniforge).
2. Download the appropriate miniforge installer for your platform from the [miniforge releases page](https://github.com/conda-forge/miniforge/releases). For example, on a Linux x86_64 machine:

`wget https://github.com/conda-forge/miniforge/releases/download/24.9.2-0/Miniforge3-24.9.2-0-Linux-x86_64.sh`

For Mac with Apple Silicon (M1/M2/M3/M4):

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
pharokka install -h
pharokka -h
```

