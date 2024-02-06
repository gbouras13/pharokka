import collections
import os

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

# from itertools import chain


def run_pyhmmer(db_dir, out_dir, threads, gene_predictor, evalue):
    """
    Runs PyHMMER on phrogs
    :param db_dir: database path
    :param out_dir: output directory
    :params threads: threads
    :param gene_predictor: phanotate or prodigal
    :param evalue: evalue threshold for pyhmmer
    :return: best_results - dictionary of top HMM hits for each protein
    """

    # define the amino acid FASTA
    amino_acid_fasta_name = gene_predictor + "_aas_tmp.fasta"
    amino_acid_fasta_file = os.path.join(out_dir, amino_acid_fasta_name)

    # define result
    Result = collections.namedtuple(
        "Result", ["protein", "phrog", "bitscore", "evalue"]
    )

    # run hmmscan and get all results
    results = []
    with pyhmmer.plan7.HMMFile(os.path.join(db_dir, "all_phrogs.h3m")) as hmms:  # hmms
        with pyhmmer.easel.SequenceFile(
            amino_acid_fasta_file, digital=True
        ) as seqs:  # amino acid sequences
            for hits in pyhmmer.hmmer.hmmscan(
                seqs, hmms, cpus=int(threads), E=float(evalue)
            ):  # run hmmscan
                protein = hits.query_name.decode()  # get protein from the hit
                for hit in hits:
                    if hit.included:
                        # include the hit to the result collection
                        results.append(
                            Result(protein, hit.name.decode(), hit.score, hit.evalue)
                        )

    # get  best results for each protein
    best_results = {}
    keep_protein = set()
    for result in results:
        if result.protein in best_results:
            previous_bitscore = best_results[result.protein].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.protein] = result
                keep_protein.add(result.protein)
            elif result.bitscore == previous_bitscore:
                if best_results[result.protein].phrog != hit.name.decode():
                    keep_protein.remove(result.protein)
        else:
            best_results[result.protein] = result
            keep_protein.add(result.protein)

    return best_results
