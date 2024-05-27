import argparse
import os
import shutil
import subprocess as sp
from argparse import RawTextHelpFormatter

import pyrodigal
import pyrodigal_gv
from Bio import SeqIO
from loguru import logger
from util import get_version


def get_input():
    parser = argparse.ArgumentParser(
        description="pharokka: fast phage annotation program",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--infile", action="store", help="Input genome file in fasta format."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Directory to write the output to.",
        default=os.path.join(os.getcwd(), "output/"),
    )
    parser.add_argument(
        "-d",
        "--database",
        action="store",
        help="Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.",
        default="Default",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads. Defaults to 1.",
        action="store",
        default=str(1),
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        action="store",
        help="Prefix for output files. This is not required.",
        default="Default",
    )
    parser.add_argument(
        "-l",
        "--locustag",
        action="store",
        help="User specified locus tag for the gff/gbk files. This is not required. A random locus tag will be generated instead.",
        default="Default",
    )
    parser.add_argument(
        "-g",
        "--gene_predictor",
        action="store",
        help='User specified gene predictor. Use "-g phanotate" or "-g prodigal" or "-g prodigal-gv" or "-g genbank". \nDefaults to phanotate usually and prodigal-gv in meta mode.',
        default="default",
    )
    parser.add_argument(
        "-m",
        "--meta",
        help="meta mode for metavirome input samples",
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--split",
        help="split mode for metavirome samples. -m must also be specified. \nWill output separate split FASTA, gff and genbank files for each input contig.",
        action="store_true",
    )
    parser.add_argument(
        "-c",
        "--coding_table",
        help="translation table for prodigal. Defaults to 11.",
        action="store",
        default="11",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        help="E-value threshold for MMseqs2 database PHROGs, VFDB and CARD and PyHMMER PHROGs database search. Defaults to 1E-05.",
        action="store",
        default="1E-05",
    )
    parser.add_argument(
        "--fast",
        "--hmm_only",
        help="Runs PyHMMER (HMMs) with PHROGs only, not MMseqs2 with PHROGs, CARD or VFDB. \nDesigned for phage isolates, will not likely be faster for large metagenomes.",
        action="store_true",
    )
    parser.add_argument(
        "--mmseqs2_only",
        help="Runs MMseqs2 with PHROGs, CARD and VFDB only (same as Pharokka v1.3.2 and prior). Default in meta mode.",
        action="store_true",
    )
    parser.add_argument(
        "--meta_hmm",
        help="Overrides --mmseqs2_only in meta mode. Will run both MMseqs2 and PyHMMER.",
        action="store_true",
    )
    parser.add_argument(
        "--dnaapler",
        help="Runs dnaapler to automatically re-orient all contigs to begin with terminase large subunit if found. \nRecommended over using '--terminase'.",
        action="store_true",
    )
    parser.add_argument(
        "--custom_hmm",
        help="Run pharokka with a custom HMM profile database suffixed .h3m. \nPlease use create this with the create_custom_hmm.py script.",
        action="store",
        default="",
    )
    parser.add_argument(
        "--genbank",
        help="Flag denoting that -i/--input is a genbank file instead of the usual FASTA file. \n The CDS calls in this file will be preserved and re-annotated.",
        action="store_true",
    )
    parser.add_argument(
        "--terminase",
        help="Runs terminase large subunit re-orientation mode. \nSingle genome input only and requires --terminase_strand and --terminase_start to be specified.",
        action="store_true",
    )
    parser.add_argument(
        "--terminase_strand",
        help='Strand of terminase large subunit. Must be "pos" or "neg".',
        action="store",
        default="nothing",
    )
    parser.add_argument(
        "--terminase_start",
        help="Start coordinate of the terminase large subunit.",
        action="store",
        default="nothing",
    )
    parser.add_argument(
        "--skip_extra_annotations",
        help="Skips tRNAscan-SE 2, MinCED and Aragorn.",
        action="store_true",
    ),
    parser.add_argument(
        "--skip_mash",
        help="Skips running mash to find the closest match for each contig in INPHARED.",
        action="store_true",
    )
    parser.add_argument(
        "--minced_args",
        help='extra commands to pass to MINced (please omit the leading hyphen for the first argument). You will need to use quotation marks e.g. --minced_args "minNR 2 -minRL 21"',
        default="",
        type=str,
    )
    parser.add_argument(
        "--mash_distance",
        help="mash distance for the search against INPHARED. Defaults to 0.2.",
        default=0.2,
        type=float,
    )
    parser.add_argument(
        "-V",
        "--version",
        help="Print pharokka Version",
        action="version",
        version=get_version(),
    )
    parser.add_argument(
        "--citation", help="Print pharokka Citation", action="store_true"
    )
    args = parser.parse_args()

    return args


def instantiate_dirs(output_dir, meta, force):
    # remove outdir on force
    if force == True:
        if os.path.isdir(output_dir) == True:
            logger.info(
                f"Removing output directory {output_dir} as -f or --force was specified."
            )
            shutil.rmtree(output_dir)

        elif os.path.isfile(output_dir) == True:
            os.remove(output_dir)
        else:
            logger.info(
                f"--force was specified even though the output directory {output_dir} does not already exist. Continuing."
            )
    else:
        if os.path.isdir(output_dir) == True or os.path.isfile(output_dir) == True:
            logger.error(
                f"The output directory {output_dir} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )

    # instantiate outdir
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
    mmseqs_dir = os.path.join(output_dir, "mmseqs/")
    if os.path.isdir(mmseqs_dir) == False:
        os.mkdir(mmseqs_dir)
    vfdb_dir = os.path.join(output_dir, "VFDB/")
    if os.path.isdir(vfdb_dir) == False:
        os.mkdir(vfdb_dir)
    CARD_dir = os.path.join(output_dir, "CARD/")
    if os.path.isdir(CARD_dir) == False:
        os.mkdir(CARD_dir)

    # tmp dir for  meta mode trnascan and phanotate
    input_tmp_dir = os.path.join(output_dir, "input_split_tmp/")
    if meta == True:
        if os.path.isdir(input_tmp_dir) == False:
            os.mkdir(input_tmp_dir)

    return output_dir


def validate_fasta(filename):
    if os.path.isfile(filename) == False:  # if file doesnt exist
        logger.error(
            f"Error: Input file {filename} does not exist. Please check your input."
        )
    else:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            logger.info(f"Checking input {filename}.")
            if any(fasta):
                logger.info(f"Input {filename} is in FASTA format.")
            else:
                logger.error("Error: Input file is not in the FASTA format.")

    # check for duplicate headers
    logger.info(f"Checking input {filename} for duplicate FASTA headers.")
    check_duplicate_headers(filename)
    logger.info(f"All headers in {filename} are unique.")


def check_duplicate_headers(fasta_file):
    """
    checks if there are duplicated in the FASTA header
    in response to Tina's issue
    https://github.com/gbouras13/pharokka/issues/293
    """
    header_set = set()

    # Iterate through the FASTA file and check for duplicate headers
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        if header in header_set:
            logger.error(
                f"Duplicate contig header found: {header}"
            )  # errors if duplicate header found
        if "#" in header:
            logger.error(
                f"# character found in contig: {header} - please remove all '#' from your contig headers"
            )  # errors if duplicate header found   
        else:
            header_set.add(header)
    # if it finished it will be fine


def validate_gene_predictor(gene_predictor, genbank_flag):
    if gene_predictor == "phanotate":
        logger.info("Phanotate will be used for gene prediction.")
    elif gene_predictor == "prodigal":
        logger.info(
            "Prodigal implemented with pyrodigal will be used for gene prediction."
        )
    elif gene_predictor == "prodigal-gv":
        logger.info(
            "Prodigal-gv implemented with pyrodigal-gv will be used for gene prediction."
        )
    elif gene_predictor == "genbank":
        if genbank_flag is False:
            logger.error(
                "Error: you have specified -g genbank without --genbank. Please just use --genbank."
            )
    else:
        logger.error(
            "Error: gene predictor was incorrectly specified. Please use 'phanotate', 'prodigal' or 'prodigal-gv'."
        )


def validate_meta(filepath_in, meta, split, genbank):
    if genbank is False:
        num_fastas = len([1 for line in open(filepath_in) if line.startswith(">")])
        if meta == True:
            if num_fastas < 2:
                logger.error(
                    "ERROR: -m meta mode specified when the input file only contains 1 contig. Please re-run without specifying -m."
                )
            else:
                message = f"{num_fastas} input contigs detected."
                logger.info(message)
                if split == True:
                    message = "Split mode activtated. Separate output FASTA, gff and genbank files will be output for each contig."
                    logger.info(message)
        else:
            if num_fastas > 1:
                message = "More than one contig detected in the input file. Re-running pharokka with -m meta mode is recommended unless this is a fragmented isolate genome. Continuing."
                logger.info(message)
            if split == True:
                message = "-s or --split was specified without -m or --meta and will be ignored. Please specify -s with -m if you want to run split mode. Continuing."
                logger.info(message)


def validate_strand(strand):
    if strand != "pos" and strand != "neg":
        logger.error(
            "Error: terminase strand was incorrectly specified. It should be either 'pos' or 'neg'. Please check your input and try again. \n"
        )


def validate_terminase_start(terminase_start):
    try:
        int(terminase_start)
    except:
        logger.error(
            "Error: terminase start coordinate specified is not an integer. Please check your input and try again. \n"
        )


def validate_terminase(filepath_in, terminase_strand, terminase_start):
    if terminase_strand == "nothing":
        logger.error(
            "Error: you specified --terminase to reorient your phage to begin with the terminase large subunit, but didn't specify its strand with --terminase_strand. Please check your input and try again. \n"
        )
    if terminase_start == "nothing":
        logger.error(
            "Error: you specified --terminase to reorient your phage to begin with the terminase large subunit, but didn't specify its start coordinate with --terminase_start. Please check your input and try again. \n"
        )
    validate_strand(terminase_strand)
    validate_terminase_start(terminase_start)
    num_fastas = len([1 for line in open(filepath_in) if line.startswith(">")])
    if num_fastas > 1:
        logger.error(
            "Error: To reorient your phage genome to begin with the terminase large subunit, you can only input 1 phage genome contig. Multiple contigs were detected. Please try again. \n"
        )


def validate_threads(threads):
    try:
        x = int(threads)
        x += 1
    except:
        message = (
            "Error: you specified a non-integer value for threads of "
            + threads
            + ". Please check your input and try Pharokka again."
        )
        logger.error(message)


#######
# dependencies
#######


def check_dependencies(skip_mash):
    """Checks the dependencies and versions
    skip_mash flag from args, won't check mash is skip mash specified
    :return:
    """
    #############
    # phanotate
    #############
    try:
        process = sp.Popen(
            ["phanotate.py", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT
        )
    except:
        logger.error("Phanotate not found. Please reinstall pharokka.")
    phan_out, _ = process.communicate()
    phanotate_version = phan_out.decode().strip()
    phanotate_major_version = int(phanotate_version.split(".")[0])
    phanotate_minor_version = int(phanotate_version.split(".")[1])
    phanotate_minorest_version = phanotate_version.split(".")[2]

    logger.info(
        f"Phanotate version found is v{phanotate_major_version}.{phanotate_minor_version}.{phanotate_minorest_version}"
    )

    if phanotate_major_version < 1:
        logger.error("Phanotate is too old - please reinstall pharokka.")
    if phanotate_minor_version < 5:
        logger.error("Phanotate is too old - please reinstall pharokka.")

    logger.info("Phanotate version is ok.")

    #############
    # mmseqs
    #############
    try:
        process = sp.Popen(["mmseqs"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("MMseqs2 not found. Please reinstall pharokka.")
    mmseqs_out, _ = process.communicate()
    mmseqs_out = mmseqs_out.decode()

    version_line = []

    for line in mmseqs_out.split("\n"):
        if "Version" in line:
            version_line.append(line)

    mmseqs_version = version_line[0].split(" ")[2]
    mmseqs_major_version = int(mmseqs_version.split(".")[0])
    mmseqs_minor_version = int(mmseqs_version.split(".")[1])

    logger.info(
        f"MMseqs2 version found is v{mmseqs_major_version}.{mmseqs_minor_version}"
    )

    if mmseqs_major_version != 13:
        logger.error("MMseqs2 is the wrong version. Please install v13.45111")
    if mmseqs_minor_version != 45111:
        logger.error("MMseqs2 is the wrong version. Please install v13.45111")

    logger.info("MMseqs2 version is ok.")

    #############
    # trnascan
    #############
    try:
        process = sp.Popen(["tRNAscan-SE", "-h"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("tRNAscan-SE not found. Please reinstall pharokka.")

    trna_out, _ = process.communicate()
    trna_out = trna_out.decode()

    version_line = []

    for line in trna_out.split("\n"):
        if "tRNAscan-SE" in line:
            version_line.append(line)

    trna_version = version_line[0].split(" ")[1]
    trna_major_version = int(trna_version.split(".")[0])
    trna_minor_version = int(trna_version.split(".")[1])
    trna_minorest_version = int(trna_version.split(".")[2])

    logger.info(
        f"tRNAscan-SE version found is v{trna_major_version}.{trna_minor_version}.{trna_minorest_version}"
    )

    if trna_major_version != 2:
        logger.error("tRNAscan-SE is the wrong version. Please re-install pharokka.")
    if trna_minor_version != 0:
        logger.error("tRNAscan-SE is the wrong version. Please re-install pharokka.")
    if trna_minorest_version < 9:
        logger.error("tRNAscan-SE is the wrong version. Please re-install pharokka.")

    logger.info("tRNAscan-SE version is ok.")

    #############
    # minced
    #############
    try:
        process = sp.Popen(["minced", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("MinCED not found. Please reinstall pharokka.")

    minced_out, _ = process.communicate()
    minced_out = minced_out.decode()

    version_line = []

    for line in minced_out.split("\n"):
        if "minced" in line:
            version_line.append(line)

    minced_version = version_line[0].split(" ")[1]
    minced_major_version = int(minced_version.split(".")[0])
    minced_minor_version = int(minced_version.split(".")[1])
    minced_minorest_version = int(minced_version.split(".")[2])

    logger.info(
        f"MinCED version found is v{minced_major_version}.{minced_minor_version}.{minced_minorest_version}"
    )

    if minced_major_version != 0:
        logger.error("MinCED is the wrong version. Please re-install pharokka.")
    if minced_minor_version != 4:
        logger.error("MinCED is the wrong version. Please re-install pharokka.")
    if minced_minorest_version < 2:
        logger.error("MinCED is the wrong version. Please re-install pharokka.")

    logger.info("MinCED version is ok.")

    #############
    # aragorn
    #############
    try:
        process = sp.Popen(["aragorn", "-h"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("ARAGORN not found. Please reinstall pharokka.")

    aragorn_out, _ = process.communicate()
    aragorn_out = aragorn_out.decode()

    version_line = []

    for line in aragorn_out.split("\n"):
        if "Dean Laslett" in line:
            version_line.append(line)

    aragorn_version = version_line[0].split(" ")[1]
    # strip off v
    aragorn_version = aragorn_version[1:]
    aragorn_major_version = int(aragorn_version.split(".")[0])
    aragorn_minor_version = int(aragorn_version.split(".")[1])
    aragorn_minorest_version = int(aragorn_version.split(".")[2])

    logger.info(
        f"ARAGORN version found is v{aragorn_major_version}.{aragorn_minor_version}.{aragorn_minorest_version}"
    )

    if aragorn_major_version != 1:
        logger.error("ARAGORN is the wrong version. Please re-install pharokka.")
    if aragorn_minor_version != 2:
        logger.error("ARAGORN is the wrong version. Please re-install pharokka.")
    if aragorn_minorest_version < 41:
        logger.error("ARAGORN is the wrong version. Please re-install pharokka.")

    logger.info("ARAGORN version is ok.")

    #############
    # mash
    #############
    if skip_mash is False:
        try:
            process = sp.Popen(["mash", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT)
        except:
            logger.error("mash not found. Please reinstall pharokka.")

        mash_out, _ = process.communicate()
        mash_out = mash_out.decode().strip()

        mash_major_version = int(mash_out.split(".")[0])
        mash_minor_version = int(mash_out.split(".")[1])

        logger.info(f"mash version found is v{mash_major_version}.{mash_minor_version}")

        if mash_major_version != 2:
            logger.error("mash is the wrong version. Please re-install pharokka.")
        if mash_minor_version < 2:
            logger.error("mash is the wrong version. Please re-install pharokka.")

        logger.info("mash version is ok.")

    #############
    # dnaapler
    #############
    try:
        process = sp.Popen(["dnaapler", "--version"], stdout=sp.PIPE, stderr=sp.STDOUT)
    except:
        logger.error("Dnaapler not found. Please reinstall pharokka.")

    dnaapler_out, _ = process.communicate()
    dnaapler_out = dnaapler_out.decode()

    dnaapler_out = dnaapler_out.split("n ")[1]

    dnaapler_major_version = int(dnaapler_out.split(".")[0])
    dnaapler_minor_version = int(dnaapler_out.split(".")[1])
    dnaapler_minorest_version = int(dnaapler_out.split(".")[2])

    logger.info(
        f"Dnaapler version found is v{dnaapler_major_version}.{dnaapler_minor_version}.{dnaapler_minorest_version}"
    )

    if dnaapler_minor_version < 2:
        logger.error("Dnaapler is the wrong version. Please re-install pharokka.")

    logger.info("Dnaapler version is ok.")

    #######
    # pyrodigal
    #######

    pyrodigal_version = pyrodigal.__version__
    pyrodigal_major_version = int(pyrodigal_version.split(".")[0])

    if pyrodigal_major_version < 3:
        logger.error("Pyrodigal is the wrong version. Please re-install pharokka.")

    logger.info(f"Pyrodigal version is v{pyrodigal_version}")
    logger.info(f"Pyrodigal version is ok.")

    #######
    # pyrodigal gv
    #######

    pyrodigal_gv_version = pyrodigal_gv.__version__
    pyrodigal_major_version = int(pyrodigal_gv_version.split(".")[0])

    if pyrodigal_major_version < 0:
        logger.error("Pyrodigal_gv is the wrong version. Please re-install pharokka.")

    logger.info(f"Pyrodigal_gv version is v{pyrodigal_gv_version}")
    logger.info(f"Pyrodigal_gv version is ok.")

    return (
        phanotate_version,
        pyrodigal_version,
        pyrodigal_gv_version,
        trna_version,
        aragorn_version,
        minced_version,
    )


def instantiate_split_output(out_dir, split):
    """
    if the split flag is true, will create these output directories
    """
    if split == True:
        single_gff_dir = os.path.join(out_dir, "single_gffs")
        single_gbk_dir = os.path.join(out_dir, "single_gbks")

        if os.path.isdir(single_gff_dir) == False:
            os.mkdir(single_gff_dir)
        if os.path.isdir(single_gbk_dir) == False:
            os.mkdir(single_gbk_dir)


def validate_custom_hmm(filename):
    suffix = ".h3m"

    logger.info(f"Checking custom hmm profile {filename}")

    if filename.endswith(suffix) is True:
        logger.info(f"{filename} checked.")
    else:
        logger.exit(
            f"{filename} does not end with .h3m . Please check your --custom_hmm parameter or use create_custom_hmm.py to create a custom HMM profile."
        )


def validate_and_extract_genbank(filename, out_dir):
    # Open the GenBank file
    try:
        with open(filename, "r") as genbank_file:
            for record in SeqIO.parse(genbank_file, "genbank"):
                output_fasta_path = f"{out_dir}/genbank.fasta"

                # Write the sequence in FASTA format to the output file
                with open(output_fasta_path, "a") as output_fasta:
                    SeqIO.write(record, output_fasta, "fasta")
    except:
        logger.error(f"{filename} is not a genbank format file.")
