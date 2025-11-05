#!/usr/bin/env python3

"""Script to convert pharokka genbank file to a form with functions parsable with Snapgene (adds the product to the ID)

python3 convert_genbank.py -i pharokka.gbk -o snapgene.gbk

"""

import argparse
import os
import sys
from argparse import RawTextHelpFormatter

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from loguru import logger


def get_input():
    """gets input for create_custom_hmm.py
    :return: args
    """
    parser = argparse.ArgumentParser(
        description="convert_genbank.py: Script to convert pharokka genbank file to a form with functions parsable with Snapgene.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        help="Input genbank file (direct from pharokka.gbk).",
    )

    parser.add_argument(
        "-o",
        "--output",
        action="store",
        default="",
        help="Input genbank file for use in snapgene.",
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    args = parser.parse_args()

    return args


# Define a function to update the record IDs
def update_records(records):
    updated_records = []
    for record in records:
        # Create an empty list to store the new features
        new_features = []

        for feature in record.features:
            if feature.type == "CDS":
                # Access specific qualifiers if needed
                id = feature.qualifiers.get("ID", ["N/A"])[0]
                product = feature.qualifiers.get("product", ["N/A"])[0]

                new_id = f"{id}_{product}"
                feature.qualifiers["id"] = [new_id]

                # Create a new feature with the updated 'gene' qualifier
                new_feature = SeqFeature(
                    location=FeatureLocation(
                        feature.location.start, feature.location.end
                    ),
                    type=feature.type,
                    qualifiers={
                        "ID": [new_id],
                        "function": feature.qualifiers.get("function"),
                        "locus_tag": [new_id],
                        "frame": feature.qualifiers.get("frame"),
                        "phrog": feature.qualifiers.get("phrog"),
                        "product": feature.qualifiers.get("product"),
                        "score": feature.qualifiers.get("score"),
                        "source": feature.qualifiers.get("source"),
                        "top_hit": feature.qualifiers.get("top_hit"),
                        "transl_table": feature.qualifiers.get("transl_table"),
                        "translation": feature.qualifiers.get("translation"),
                    },
                )
                new_features.append(new_feature)
            else:
                new_features.append(feature)

        # Add the new feature to the list of new features

        record.features = new_features
        updated_records.append(record)
    return updated_records


def main():
    args = get_input()

    logger.add(lambda _: sys.exit(1), level="ERROR")

    #### force
    if args.force == True:
        if os.path.isfile(args.output) == True:
            logger.info(
                f"Removing output file {args.output} as -f or --force was specified."
            )
            os.remove(args.output)
        else:
            logger.info(
                f"--force was specified even though the output file {args.output} does not already exist. Continuing."
            )
    else:
        if os.path.isfile(args.output) == True:
            logger.error(
                f"The output file {args.output} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )

    logger.info(
        f"Creating a converted Genbank file for use with Snapgene {args.output}."
    )

    # Read the GenBank file and parse the records
    with open(args.input.strip(), "rt") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()

    # Update the record IDs
    updated_records = update_records(records)

    # Write the updated records to a new GenBank file
    SeqIO.write(updated_records, args.output, "genbank")


if __name__ == "__main__":
    main()
