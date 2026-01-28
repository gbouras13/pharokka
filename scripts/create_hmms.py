# coding: utf-8
"""Script to create hmm profiles for each PHROG with pyhmmer

Note:
    Requires all MSA (from the PHROG website) to be in the working directory in the 'MSA_Phrogs_M50_FASTA' directory.

python3 create_hmms.py

"""

import os
from pathlib import Path

import pyhmmer

alphabet = pyhmmer.easel.Alphabet.amino()
background = pyhmmer.plan7.Background(alphabet)


MSA_dir = "MSA_Phrogs_M50_FASTA"

HMM_dir = "MSA_Phrogs_M50_HMM"

# Check if the directory already exists
if not os.path.exists(HMM_dir):
    # Create the output directory
    os.mkdir(HMM_dir)

# for when they update PHROGs
number_of_phrogs = 38880

# loop over each PHROG
for i in range(1, number_of_phrogs + 1):
    # read in each msa
    with pyhmmer.easel.MSAFile(
        f"{MSA_dir}/phrog_{i}.fma", digital=True, alphabet=alphabet
    ) as msa_file:
        msa = msa_file.read()

    name = f"phrog_{i}"
    print(name)
    # convert to bytes
    msa.name = name

    # build the MSA
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)

    # writes all the PHROGS
    with open(f"{HMM_dir}/phrog_{i}.hmm", "wb") as output_file:
        hmm.write(output_file)


# to concatenate all hmms

hmms = []

# Specify the directory path
HMM_dir = Path("MSA_Phrogs_M50_HMM")

# reads and saves the hmms
for i in range(1, 38881):
    f = f"{HMM_dir}/phrog_{i}.hmm"
    with pyhmmer.plan7.HMMFile(f) as hmm_file:
        hmm = hmm_file.read()
    hmms.append(hmm)
    # to let you know where you are at
    print(i)

# writes all out together to .h3m, .h3p, .h3i, .h3f files prefixed "all_phrogs"
pyhmmer.hmmer.hmmpress(hmms, "all_phrogs")
