import os

from setuptools import setup


# recursively load package files
def package_files(directory):
    paths = []
    for path, _, filenames in os.walk(directory):
        for filename in filenames:
            if not filename.endswith(".py"):
                paths.append(os.path.join("..", path, filename))
    return paths


# read long description
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="Pharokka",
    version="1.7.2",
    author="George Bouras",
    author_email="george.bouras@adelaide.edu.au",
    description="Fast phage annotation tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gbouras13/pharokka",
    scripts=[
        "bin/pharokka.py",
        "bin/pharokka_multiplotter.py",
        "bin/pharokka_plotter.py",
        "bin/install_databases.py",
        "bin/pharokka_proteins.py",
        "bin/citation.py",
        "bin/databases.py",
        "bin/external_tools.py",
        "bin/input_commands.py",
        "bin/hmm.py",
        "bin/plot.py",
        "bin/post_processing.py",
        "bin/processes.py",
        "bin/proteins.py",
        "bin/util.py",
        "bin/version.py",
        "bin/create_custom_hmm.py",
        "bin/custom_db.py",
    ],
    packages=["pharokka_runner"],
    package_dir=dict(pharokka_runner="bin"),
    package_data=dict(pharokka_runner=package_files("bin/")),
    include_package_data=True,
    license="MIT License",
    platforms=["Unix"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.5",
    install_requires=[
        "setuptools>=67.7.2",
        "loguru>=0.5.4",
        "pyyaml>=6.0",
        "pandas>=1.4.2",
        "biopython>=1.80",
        "pyhmmer>=0.10.0",
        "black>=22.3.0",
        "isort>=5.10.1",
        "pytest>=6.2.5",
        "pytest-cov>=3.0.0",
        "alive-progress>=3.0.1",
        "requests>=2.25.1",
        "bcbio-gff>=0.7.0",
        "pyrodigal>=3.0.0",
        "pyrodigal_gv>=0.1.0",
    ],
)
