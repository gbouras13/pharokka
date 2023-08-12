import contextlib
import itertools
import os
import typing

import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile, HMM

from pyhmmer.plan7 import HMMFile


# from itertools import chain

import os
from pathlib import Path




# # class HMMFiles():
# #     def __init__(self, files):
# #         self.hmmfiles = list(files)

#     # def __enter__(self):
#     #     return chain.from_iterable(self.hmmfiles)

#     # def __exit__(self, *args):
#     #     for f in self.hmmfiles:
#     #         f.close()

# # class HMMFiles(typing.Iterable[HMM]):
# #     def __init__(self, files):
# #         self.files = list(files)
# #     def __iter__(self):
# #         for file in self.files:
# #             with HMMFile(file) as hmm_file:
# #                 yield from hmm_file
# #     def __enter__(self) -> typing.Iterable[HMM]:
# #         return chain.from_iterable(self.files)
# #     def __exit__(self, *args):
# #         for f in self.files:
# #             f.close()





# # class HMMFiles(typing.ContextManager[typing.Iterable[HMM]]):
# #     def __init__(self, files: list['os.PathLike[bytes]']) -> None:
# #         self.stack = contextlib.ExitStack()
# #         self.hmmfiles = [self.stack.enter_context(HMMFile(f)) for f in files]

# #     def __enter__(self) -> typing.Iterable[HMM]:
# #         return itertools.chain.from_iterable(self.hmmfiles)

# #     def __exit__(self, exc_value: object, exc_type: object, traceback: object) -> None:
# #         self.stack.close()





# # with SequenceFile("MZ130495/prodigal.faa", digital=True) as sequences:
# #     targets = sequences.read_block()






# # Specify the suffix you're looking for
# suffix = '.hmm'



# # List files with the specified suffix
# files_with_suffix = [file for file in HMM_dir.iterdir() ]

# from pathlib import Path
# import glob

# directory_path = Path("/path/to/your/directory")
# pattern = "*.hmm"  # Replace with your desired file pattern

# files = HMM_dir.glob(pattern)

# # files = list(HMM_dir.glob(pattern))

# directory_path = "MSA_Phrogs_M50_HMM"
# files_in_directory = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]

# print(files)

# class HMMFiles(typing.Iterable[HMM]):
#     def __init__(self, files):
#         self.files = list(files)
#     def __iter__(self):
#         for file in self.files:
#             with HMMFile(file) as hmm_file:
#                 yield from hmm_file
#     def __enter__(self) -> typing.Iterable[HMM]:
#         return chain.from_iterable(self.files)
#     def __exit__(self, *args):
#         for f in self.files:
#             f.close()




import time
import pyhmmer

# with pyhmmer.plan7.HMMFile("all_phrogs.h3m") as hmms:
#     with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True) as seqs:
#         t1 = time.time()
#         total = sum(len(hits) for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs))
#         print(f"- hmmsearch found a total of {total} hits in {time.time() - t1:.3} seconds")

with pyhmmer.plan7.HMMFile("all_phrogs.h3m") as hmms:
    with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True) as seqs:
        t1 = time.time()
        total = sum(len(hits) for hits in pyhmmer.hmmer.hmmscan(seqs, hmms, E=1e-02))
        print(f"- hmmscan found a total of {total} hits in {time.time() - t1:.3} seconds")


import collections
Result = collections.namedtuple("Result", ["query", "cog", "bitscore", "evalue"])


hits = pyhmmer.hmmer.hmmsearch(hmms, seqs)


import collections
Result = collections.namedtuple("Result", ["protein", "phrog", "bitscore", "evalue"])

with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True) as seqs_file:
    proteins = seqs_file.read_block()

# cpus = threads
# get all results
results = []
with pyhmmer.plan7.HMMFile("all_phrogs.h3m") as hmms:
    with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True) as seqs:
        for hits in pyhmmer.hmmer.hmmscan(seqs, hmms, cpus=5, E=1e-02):
            phrog = hits.query_name.decode()
            for hit in hits:
                if hit.included:
                    # include the 
                    results.append(Result( phrog, hit.name.decode(), hit.score, hit.evalue))




# print(results)

# get  best restuls
best_results = {}
keep_protein = set()
for result in results:
    if result.protein in best_results:
        previous_bitscore = best_results[result.protein].bitscore
        if result.bitscore > previous_bitscore:
            best_results[result.protein] = result
            keep_protein.add(result.protein)
        elif result.bitscore == previous_bitscore:
            if best_results[result.protein].phrog != hit.phrog:
                keep_protein.remove(result.protein)
    else:
        best_results[result.protein] = result
        keep_protein.add(result.protein)


print(best_results)


# import pyhmmer.easel
# with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True) as seqs_file:
#     proteins = seqs_file.read_block()

# import collections
# Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])


# results = []
# for hits in pyhmmer.hmmsearch(hmms, proteins, bit_cutoffs="trusted"):
#     cog = hits.query_name.decode()
#     for hit in hits:
#         if hit.included:
#             results.append(Result(hit.name.decode(), cog, hit.score))
# print(results)


# class HMMFiles(typing.Iterable[HMM]):
#     def __init__(self, files):
#         self.files = list(files)
#     def __iter__(self):
#         for file in self.files:
#             with HMMFile(file) as hmm_file:
#                 yield from hmm_file
#     def __enter__(self):
#         return chain.from_iterable(self.files)
#     def __exit__(self, *args):
#         for f in self.files:
#             f.close()



# with SequenceFile("MZ130495/prodigal.faa", digital=True) as sequences:
#     targets = sequences.read_block()

# with HMMFiles(files) as hmm_files:
#     all_hits = list(pyhmmer.hmmsearch(hmm_files, targets))




# # with HMMFiles('MSA_Phrogs_M50_HMM/phrog_16490.hmm', 'MSA_Phrogs_M50_HMM/phrog_19368.hmm') as hmm_files:
# #     all_hits = list(pyhmmer.hmmsearch(hmm_files, targets))


# print("HMMs searched:", len(all_hits))
# print("Hits found:   ", sum(len(hits) for hits in all_hits))



# # pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
# # with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True, alphabet=alphabet) as seq_file:
# #     hits = pipeline.search_hmm(hmm, seq_file)



# # with pyhmmer.plan7.HMMFile("all_phrogs.hmm") as hmm_file:
# #     hmm = hmm_file.read()

# # print(hmm)

# # with pyhmmer.easel.SequenceFile("MZ130495/prodigal.faa", digital=True) as seq_file:
# #     sequences = seq_file.read_block()

# # pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
# # hits = pipeline.search_hmm(hmm, sequences)

# # # ali = hits[0].domains[0].alignment
# # # print(ali)