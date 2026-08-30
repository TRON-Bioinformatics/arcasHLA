#!/bin/env python3
"""
Central arg-parsing and script dispatch.
"""

import argparse
import sys

from arcas_utilities import tool_is_available

from convert import build_arg_parser as build_convert_parser
from customize import build_arg_parser as build_cusomize_parser
from extract import build_arg_parser as build_extract_parser
from genotype import build_arg_parser as build_genotype_parser
from merge import build_arg_parser as build_merge_parser
from partial import build_arg_parser as build_partial_parser
from quant import build_arg_parser as build_quant_parser
from reference import build_arg_parser as build_reference_parser


def check_tool_depencencies():
    """
    Docstring for check_tool_depencencies
    """
    tools = [
        "awk",
        "bedtools",
        "cat",
        "git",
        "kallisto",
        "mkdir",
        "pigz",
        "rm",
        "samtools",
        "zcat",
    ]

    missing_tools = []

    for tool in tools:
        if not tool_is_available(tool):
            missing_tools.append(tool)

    if len(missing_tools) > 0:
        missing_tool_str = ", ".join(missing_tools)

        raise Exception(
            f"Can't find tools/commands {missing_tool_str}. Ensure the packages they "
            "belong to are installed. Within this repo, there *should* be a conda env "
            "file with all dependencies required to use this package."
        )


def build_arg_parser():
    parser = argparse.ArgumentParser(prog="arcasHLA")

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        description="available subcommands",
        help=(
            "extract: extracts chromosome 6 reads from bam, "
            "genotype: types HLA genes from extracted reads, "
            "partial: types partial HLA genes from extracted reads, "
            "customize: create custom HLA reference, "
            "quant: allele specific HLA quantification, "
            "merge: processes results into a tab-separated table, "
            "convert: converts HLA nomenclature/resolution, "
            "reference: build an external, versioned HLA reference."
        ),
    )
    build_convert_parser(subparsers)
    build_cusomize_parser(subparsers)
    build_extract_parser(subparsers)
    build_genotype_parser(subparsers)
    build_merge_parser(subparsers)
    build_partial_parser(subparsers)
    build_quant_parser(subparsers)
    build_reference_parser(subparsers)

    return parser


def main(args):
    parser = build_arg_parser()

    parsed_args = parser.parse_args(args)

    check_tool_depencencies()

    # All parsers have the `run_function` default, allowing straightforward dispatch of
    # the corresponding subparser's main run function.
    parsed_args.run_function(parsed_args)


if __name__ == "__main__":
    main(sys.argv[1:])
