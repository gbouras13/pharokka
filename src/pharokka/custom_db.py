import collections
import os

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile


def run_custom_pyhmmer(custom_hmm, out_dir, threads, gene_predictor, evalue):
    """
    Runs pyhmmer on custom database
    :param custom_hmm: custom .h3m hmm file
    :param out_dir: output directory
    :params threads: threads
    :param gene_predictor: phanotate or prodigal
    :param evalue: evalue threshold for pyhmmer
    :return: best_results - dictionary of top HMM hits for each protein
    """

    amino_acid_fasta_name = gene_predictor + "_aas_tmp.fasta"
    amino_acid_fasta_file = os.path.join(out_dir, amino_acid_fasta_name)

    Result = collections.namedtuple(
        "Result", ["protein", "custom_hmm_id", "bitscore", "evalue"]
    )

    results = []
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(custom_hmm, alphabet=alphabet) as hmms:
        with pyhmmer.easel.SequenceFile(
            amino_acid_fasta_file, digital=True, alphabet=alphabet
        ) as seqs:
            for hits in pyhmmer.hmmer.hmmscan(
                seqs, hmms, cpus=int(threads), E=float(evalue)
            ):
                protein = hits.query.name
                for hit in hits:
                    if hit.included:
                        results.append(
                            Result(protein, hit.name, hit.score, hit.evalue)
                        )

    best_results = {}
    for result in results:
        if result.protein not in best_results or result.bitscore > best_results[result.protein].bitscore:
            best_results[result.protein] = result

    return best_results
