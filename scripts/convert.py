#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------
#   convert.py: changes HLA resolution and converts nomenclature to p/g-groups.
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
#   This file is part of arcasHLA.
#
#   arcasHLA is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   arcasHLA is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with arcasHLA.  If not, see <https://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------------

import sys
import os
import json
import argparse
import pandas as pd

from argparse import RawTextHelpFormatter

import config
from arcas_utilities import process_allele

# -------------------------------------------------------------------------------


def convert_allele(allele, resolution, p_group, g_group, force=False):
    """Checks nomenclature of input allele and returns converted allele."""
    i = len(allele.split(":"))

    # Input: P-group allele
    if allele[-1] == "P":
        if resolution == "g-group":
            sys.exit("[convert] Error: p-group cannot be converted " + "to g-group.")

        # Output: 1-field allele unless forced
        elif type(resolution) == int:
            if resolution > 1 and not force:
                sys.exit(
                    "[convert] Error: p-group cannot be "
                    + "converted to %.0f fields." % resolution
                )
            allele = process_allele(allele[:-1], resolution)

    # Input: G-group allele
    elif allele[-1] == "G":

        # Output: 1-field allele unless forced
        if type(resolution) == int:
            if resolution > 1 and not force:
                sys.exit(
                    "[convert] Error: g-group cannot be converted"
                    + "to %.0f fields." % resolution
                )
            allele = process_allele(allele[:-1], resolution)

        # Output: P-group allele
        elif resolution == "p-group":
            if allele[:-1] in p_group[str(i)]:
                allele = p_group[str(i)][allele[:-1]]

            elif process_allele(allele[:-1], i - 1) in p_group[str(i)]:
                allele = p_group[str(i)][process_allele(allele[:-1], i - 1)]

    # Input: ungrouped allele
    # Output: G-group allele
    elif resolution == "g-group":
        if allele in g_group[str(i)]:
            allele = g_group[str(i)][allele]
        elif allele[-1] != "N":
            allele = process_allele(allele, 3)

    # Input: ungrouped allele
    # Output: P-group allele
    elif resolution == "p-group":
        if allele in p_group[str(i)]:
            allele = p_group[str(i)][allele]

    # Input: ungrouped allele
    # Output: reduced resolution, ungrouped allele
    elif type(resolution) == int:
        allele = process_allele(allele, resolution)

    return allele


def do_conversion(file: str, resolution_string, outfile=None, force=False):
    # p_group, g_group = pickle.load(open(hla_convert,'rb'))
    # to do, test this
    with open(config.hla_convert_json, "r") as json_file:
        p_group, g_group = json.load(json_file)

    # Check input resolution
    accepted_fields = {"1", "2", "3", "4"}
    accepted_groupings = {"g-group", "p-group"}

    resolution = None

    if resolution_string in accepted_fields:
        resolution = int(resolution_string)
    elif resolution_string.lower() in accepted_groupings:
        resolution = resolution_string.lower()

    if not resolution:
        sys.exit(
            "[convert] Error: output resolution is needed "
            + "(1, 2, 3, g-group, p-group)."
        )

    # Create outfile name
    if not outfile:
        outfile = [
            os.path.splitext(os.path.basename(file))[0],
            resolution_string.lower(),
            "tsv",
        ]
        outfile = ".".join(outfile)

    # Load input genotypes
    df_genotypes = pd.read_csv(file, sep="\t").set_index("subject")
    genotypes = df_genotypes.to_dict("index")

    for subject, genotype in genotypes.items():
        for gene, allele in genotype.items():
            if type(allele) != str:
                continue

            genotypes[subject][gene] = convert_allele(
                allele, resolution, p_group, g_group, force
            )

    pd.DataFrame(genotypes).T.rename_axis("subject").to_csv(outfile, sep="\t")


def build_arg_parser(super_parser=None, subcommand_name="convert"):
    parser_args = {
        "prog": "arcasHLA convert",
        "usage": "%(prog)s [options]",
        "add_help": False,
        "formatter_class": RawTextHelpFormatter,
    }

    if not super_parser:
        parser = argparse.ArgumentParser(**parser_args)

    else:
        parser = super_parser.add_parser(name=subcommand_name, **parser_args)

    parser.add_argument(
        "file",
        help="tsv containing HLA genotypes, see github for "
        + "example file structure.\n\n",
        type=str,
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="show this help message and exit\n\n",
        default=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-r",
        "--resolution",
        help="output resolution (1,2,3) or grouping " + "(g-group, p-group)\n\n",
        metavar="",
    )

    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="output file\n  default: " + "./file_basename.resolution.tsv\n\n",
        default="",
        metavar="",
    )

    parser.add_argument(
        "-f",
        "--force",
        help="force conversion for grouped alleles even if "
        + "it results in loss of resolution",
        action="count",
        default=False,
    )

    parser.add_argument("-v", "--verbose", action="count", default=False)

    parser.set_defaults(
        run_function=lambda parsed_args: do_conversion(
            parsed_args.file,
            parsed_args.resolution,
            parsed_args.outfile,
            parsed_args.force,
        )
    )

    return parser


# -------------------------------------------------------------------------------


def main(args):
    parser = build_arg_parser()
    parsed_args = parser.parse_args(args)

    parsed_args.run_function(parsed_args)


if __name__ == "__main__":
    main(sys.argv[1:])


# -------------------------------------------------------------------------------
