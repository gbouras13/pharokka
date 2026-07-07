import os

import pyrodigal_gv
from Bio import SeqIO


def run_pyrodiga_gv(filepath_in, out_dir, coding_table):
    """
    Gets CDS using pyrodigal_gv
    :param filepath_in: input filepath
    :param out_dir: output directory
    :param coding_table: coding table
    :return:
    """

    orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)

    with open(os.path.join(out_dir, "prodigal_out.gff"), "w") as dst:
        for i, record in enumerate(SeqIO.parse(filepath_in, "fasta")):
            genes = orf_finder.find_genes(str(record.seq))
            genes.write_gff(dst, sequence_id=record.id)

    with open(os.path.join(out_dir, "prodigal_out_tmp.fasta"), "w") as dst:
        for i, record in enumerate(SeqIO.parse(filepath_in, "fasta")):
            genes = orf_finder.find_genes(str(record.seq))
            genes.write_genes(dst, sequence_id=record.id)
