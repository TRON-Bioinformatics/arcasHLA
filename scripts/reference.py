#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------
#   reference.py: builds HLA references for genotyping.
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
import re
import json

# import pickle
import argparse
import logging as log
import numpy as np

from argparse import RawTextHelpFormatter
from os.path import isfile, isdir, dirname, realpath

from textwrap import wrap
from scipy import stats
from collections import defaultdict

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import config
from arcas_utilities import *

# -------------------------------------------------------------------------------

__version__ = "0.4.0"
__date__ = "2022-01-27"

# -------------------------------------------------------------------------------
#   Fetch and process IMGTHLA database
# -------------------------------------------------------------------------------
def get_mode(lengths):
    return stats.mode(lengths)[0]


def ensure_ref_exists():
    """Check if IMGTHLA and constructed HLA references exist."""

    if not isfile(config.hla_dat):
        force_clone_imgt_hla_db()
        create_conversion_tables(False)
        create_fasta_ref_from_hla_allele_data()


def force_clone_imgt_hla_db():
    """Clones IMGTHLA github to database."""

    if isdir(config.imgthla_dir):
        run_command(["rm", "-rf", config.imgthla_dir])

    command = ["git", "clone", config.imgthla_git, config.imgthla_dir]
    run_command(command, "[reference] Cloning IMGT/HLA database:")


def checkout_imgt_hla_db(commithash, verbose=True):
    """Checks out a specific IMGTHLA github version given a commithash."""

    if not isfile(config.hla_dat):
        force_clone_imgt_hla_db()

    command = ["git", "-C", config.imgthla_dir, "checkout", "-f", commithash]
    if verbose:
        run_command(command, "[reference] Checking out IMGT/HLA:")
    else:
        run_command(command)


def hla_dat_commit(print_version=False):
    """Returns commithash of downloaded IMGTHLA database."""

    results = run_command(["git", "-C", config.imgthla_dir, "rev-parse HEAD"])
    commit = results.stdout.decode()
    if print_version:
        log.info(commit)

    return commit[:-1]


def parse_hla_allele_data():
    """Processes IMGTHLA database, returning HLA sequences, exon locations,
    lists of complete and partial alleles and possible exon combinations.
    """

    sequences = dict()
    utrs = defaultdict(dict)
    exons = defaultdict(dict)
    gene_exons = defaultdict(set)

    sequence = partial = utr = exon = False

    gene_set = set()
    complete_alleles = set()
    complete_2fields = set()
    partial_alleles = set()

    with open(config.hla_dat, "r", encoding="UTF-8") as file:
        lines = file.read().splitlines()

    # Check if hla.dat failed to download
    if len(lines) < 10:
        sys.exit("[reference] Error: dat/IMGTHLA/hla.dat empty or corrupted.")

    for line in lines:
        # Denotes end of sequence, add allele to database
        if line.startswith("//"):
            if sequence and allele in exons:
                sequences[allele] = seq
                gene_exons[gene].add(number)
                gene_set.add(gene)

                if not partial:
                    complete_alleles.add(allele)
                    complete_2fields.add(process_allele(allele, 2))

                else:
                    partial_alleles.add(allele)
            partial = False

        # Denotes partial alleles
        elif line.startswith("FT") and "partial" in line:
            partial = True

        # Allele name and gene
        elif line.startswith("FT") and re.search('allele\="HLA-', line):
            allele = re.split("HLA-", re.sub('["\n]', "", line))[1]
            gene = get_gene(allele)

            exon = sequence = False
            seq = ""

        # Exon coordinates
        elif line.startswith("FT") and re.search("exon", line):
            info = re.split("\s+", line)
            start = int(info[2].split("..")[0]) - 1
            stop = int(info[2].split("..")[1])
            exon_coord = [start, stop]
            exon = True

        # Exon number on following line
        elif exon:
            number = re.split('"', line)[1]
            exons[allele][number] = exon_coord
            exon = False

        # UTRs
        elif line.startswith("FT") and (re.search("\sUTR\s", line)):
            info = re.split("\s+", line)
            start = int(info[2].split("..")[0]) - 1
            stop = int(info[2].split("..")[1])
            utr_coord = [start, stop]

            if allele not in exons:
                utrs[allele]["utr5"] = utr_coord
            else:
                utrs[allele]["utr3"] = utr_coord

        # Start of sequence
        elif line.startswith("SQ"):
            sequence = True

        elif sequence and line.startswith(" "):
            seq += "".join(line.split()[:-1]).upper()

    # select only 2-field partial alleles
    partial_alleles = {
        allele
        for allele in partial_alleles
        if process_allele(allele, 2) not in complete_2fields
    }

    # get most common final exon length to truncate stop-loss alleles
    final_exon_length = defaultdict(list)
    for allele in complete_alleles:
        gene = get_gene(allele)
        exon = sorted(gene_exons[gene])[-1]

        if exon not in exons[allele]:
            continue

        start, stop = exons[allele][exon]
        final_exon_length[gene].append(stop - start)

    for gene, lengths in final_exon_length.items():
        exon = sorted(gene_exons[gene])[-1]
        length = get_mode(lengths)
        final_exon_length[gene] = [exon, length]

    return (
        complete_alleles,
        partial_alleles,
        gene_set,
        sequences,
        utrs,
        exons,
        final_exon_length,
    )


def parse_hla_nomenclature(hla_nom):
    """Processes nomenclature files for arcasHLA convert."""
    allele_to_group = defaultdict(dict)

    single_alleles = set()
    grouped_alleles = set()

    for line in open(hla_nom, "r", encoding="UTF-8"):
        if line.startswith("#"):
            continue

        gene, alleles, group = line.split(";")
        alleles = [gene + allele for allele in alleles.split("/")]
        if len(group) == 1:
            single_alleles.add(alleles[0])
            continue

        group = gene + group[:-1]

        for allele in alleles:
            grouped_alleles.add(process_allele(allele, 2))
            for i in range(2, 5):
                allele_to_group[i][process_allele(allele, i)] = group

    # Alleles not included in a group
    for allele in single_alleles:
        # Checks if 2-field allele already represented by a group
        if process_allele(allele, 2) not in grouped_alleles:
            for i in range(2, 5):
                allele_to_group[i][process_allele(allele, i)] = process_allele(
                    allele, 2
                )
        else:
            allele_to_group[2][allele] = process_allele(allele, 3)
            allele_to_group[3][allele] = process_allele(allele, 3)
            allele_to_group[4][allele] = allele

    return allele_to_group


# -------------------------------------------------------------------------------
#   Saving reference files
# -------------------------------------------------------------------------------


def write_reference(sequences, info, fasta, idx, database, type):
    """Writes and idxes HLA references."""
    with open(fasta, "w") as file:
        SeqIO.write(sequences, file, "fasta")

    commithash = hla_dat_commit()
    # with open(database,'wb') as file:
    #    pickle.dump([commithash,info],file)
    with open(database, "w") as file:
        if len(info) == 4:
            json.dump(
                [
                    commithash,
                    [
                        list(info[0]),
                        json.dumps(info[1], cls=NumpyEncoder),
                        json.dumps(info[2], cls=NumpyEncoder),
                        json.dumps(info[3], cls=NumpyEncoder),
                    ],
                ],
                file,
            )
        if len(info) == 6:
            json.dump(
                [
                    commithash,
                    [
                        list(info[0]),
                        json.dumps(info[1], cls=NumpyEncoder),
                        json.dumps(info[2], cls=NumpyEncoder),
                        json.dumps(info[3], cls=NumpyEncoder),
                        json.dumps(info[4], cls=NumpyEncoder),
                        list(info[5]),
                    ],
                ],
                file,
            )

    run_command(
        ["kallisto", "index", "-i", idx, fasta],
        "[reference] indexing " + type + " reference with Kallisto:",
    )


# -------------------------------------------------------------------------------
#   Constructing reference
# -------------------------------------------------------------------------------


def get_exon_combinations():
    """Generates exon combinations used in partial allele typing."""
    exon_combinations = []
    exon_set = set()
    for i in range(2, 8):
        exon_set |= {str(i)}
        exon_combinations.append(sorted(exon_set))
        if i > 2:
            exon_combinations.append(sorted(exon_set | {"1"}))
    return exon_combinations


def create_fasta_ref_from_hla_allele_data():
    """Constructs HLA reference from processed sequences and exon locations."""

    log.info("[reference] IMGT/HLA database version:\n")
    hla_dat_commit(True)

    log.info("[reference] Processing IMGT/HLA database")

    # Constructs cDNA sequences for alleles and adds UTRs to the set of
    # non-coding sequences
    def build_complete(allele):
        gene = get_gene(allele)
        allele_exons = sorted(exons[allele].items())
        coords = [[start, stop] for n, (start, stop) in allele_exons]
        seq = [sequences[allele][start:stop] for start, stop in coords]
        seq = "".join(seq)

        exon, exon_length = final_exon_length[gene]
        if exon in exons[allele]:
            start, stop = exons[allele][exon]
            if stop - start > exon_length:
                seq = seq[: exon_length - (stop - start) + 1]

        cDNA[seq].add(allele)
        gene_length[gene].append(len(seq))

        if allele in utrs:
            for start, stop in utrs[allele].values():
                seq = sequences[allele][start:stop]
                other.add(seq)

    # Constructs exon combination sequences for complete and
    # partial alleles
    def build_combination(allele):
        for exon_group in exon_combinations:
            if set(exon_group) & set(exons[allele]) != set(exon_group):
                continue

            coords = []
            for n in exon_group:
                coords.append(exons[allele][n])
            seq = [sequences[allele][start:stop] for start, stop in coords]
            seq = "".join(seq)
            combo[str(exon_group)][seq].add(allele)

    # Adds cDNA sequences of complete alleles and separate UTRs to a
    # list of sequence records
    def complete_records(cDNA, other):
        seq_out = []
        allele_idx = dict()
        lengths = dict()

        # Adds coding sequences
        cDNA = sorted(cDNA.items(), key=lambda x: x[1])
        for i, (seq, alleles) in enumerate(cDNA):
            idx = str(i)
            record = SeqRecord(Seq(seq), id=str(idx), description="")
            seq_out.append(record)

            allele_idx[idx] = sorted(alleles)
            lengths[idx] = len(seq)

        # Adds UTRs
        offset = i + 1
        for i, seq in enumerate(other):
            idx = str(i + offset)

            record = SeqRecord(Seq(seq), id=str(idx), description="")
            seq_out.append(record)

            allele_idx[idx] = None

        return seq_out, allele_idx, lengths

    # Adds exon combination sequences of complete alleles to list of
    # sequence records, including UTRs
    def partial_records(sequences, other):
        seq_out = []
        exon_idx = dict()
        allele_idx = dict()
        lengths = dict()

        # Adds exon combination sequences
        offset = 0
        for exon in sorted(sequences):
            exon_sequences = sorted(sequences[exon].items(), key=lambda x: x[1])
            for i, (seq, alleles) in enumerate(exon_sequences):
                idx = str(i + offset)
                length = len(seq)

                record = SeqRecord(Seq(seq), id=str(idx), description="")

                seq_out.append(record)

                allele_idx[idx] = sorted(alleles)
                lengths[idx] = len(seq)
                exon_idx[idx] = exon

            offset += i + 1

        # Adds UTRs
        for i, seq in enumerate(other):
            idx = str(i + offset)

            record = SeqRecord(Seq(seq), id=str(idx), description="")

            seq_out.append(record)

            allele_idx[idx] = None

        return seq_out, allele_idx, lengths, exon_idx

    (
        complete_alleles,
        partial_alleles,
        gene_set,
        sequences,
        utrs,
        exons,
        final_exon_length,
    ) = parse_hla_allele_data()

    exon_combinations = get_exon_combinations()

    gene_length = defaultdict(list)
    cDNA = defaultdict(set)
    combo = {str(i): defaultdict(set) for i in exon_combinations}
    other = set()

    # Build sequences for each allele
    for allele in complete_alleles:
        build_complete(allele)
        build_combination(allele)

    for allele in partial_alleles:
        build_combination(allele)

    cDNA = {seq: sorted(alleles) for seq, alleles in cDNA.items()}
    other = list(other)
    gene_length = {g: get_mode(lengths) for g, lengths in gene_length.items()}

    log.info("[reference] Building HLA database")
    seq_out, allele_idx, lengths = complete_records(cDNA, other)
    write_reference(
        seq_out,
        [gene_set, allele_idx, lengths, gene_length],
        # hla_fa, hla_idx, hla_p,
        config.hla_fa,
        config.hla_idx,
        config.hla_json,
        "complete",
    )

    log.info("[reference] Building partial HLA database")
    seq_out, allele_idx, lengths, exon_idx = partial_records(combo, other)
    partial_exons = {allele: exons[allele] for allele in partial_alleles}
    write_reference(
        seq_out,
        [gene_set, allele_idx, exon_idx, lengths, partial_exons, partial_alleles],
        # partial_fa, partial_idx, partial_p,
        config.partial_fa,
        config.partial_idx,
        config.partial_json,
        "partial",
    )


def create_conversion_tables(reset=False):
    """Creates conversion tables for arcasHLA convert."""

    log.info("[reference] Building nomenclature conversion tables.")

    commit = hla_dat_commit()

    if reset:
        checkout_imgt_hla_db("origin", False)

    p_group = parse_hla_nomenclature(config.hla_nom_p)
    g_group = parse_hla_nomenclature(config.hla_nom_g)

    # with open(hla_convert, 'wb') as file:
    #    pickle.dump([p_group,g_group], file)
    # todo, test this:
    with open(config.hla_convert_json, "w") as file:
        json.dump([p_group, g_group], file)

    if reset:
        checkout_imgt_hla_db(commit, False)


class NumpyEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.int64):
            return int(o)
        return json.JSONEncoder.default(self, o)


# -------------------------------------------------------------------------------
#   Main
# -------------------------------------------------------------------------------

if __name__ == "__main__":

    # with open(parameters, 'rb') as file:
    #    _, _, versions = pickle.load(file)
    #    temp1, temp2, versions = pickle.load(file)
    # with open(parameters_json, 'w') as file:
    #    json.dump([list(temp1),list(temp2),versions],file)
    with open(config.parameters_json, "r") as file:
        _, _, versions = json.load(file)

    parser = argparse.ArgumentParser(
        prog="arcasHLA reference",
        usage="%(prog)s [options]",
        add_help=False,
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="show this help message and exit\n\n",
        default=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--update", action="count", help="update to latest IMGT/HLA version\n\n"
    )

    parser.add_argument("--rebuild", help="rebuild HLA database\n\n", action="count")

    parser.add_argument(
        "--version",
        type=str,
        help="checkout IMGT/HLA version using version\n"
        + "\n".join(wrap("options: " + ", ".join(sorted(versions.keys())), 60))
        + "\n\n",
        default=False,
        metavar="",
    )

    parser.add_argument(
        "--commit",
        type=str,
        help="checkout IMGT/HLA version using commithash\n\n",
        default=False,
        metavar="",
    )

    parser.add_argument("-v", "--verbose", action="count", default=False)

    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(level=log.DEBUG, format="%(message)s")
    else:
        handlers = [log.StreamHandler()]
        log.basicConfig(format="%(message)s")

    log.info("")
    hline()

    check_path(config.ref_dir)

    if args.update:
        log.info("[reference] Updating HLA reference")
        checkout_imgt_hla_db("origin")
        create_conversion_tables(False)
        create_fasta_ref_from_hla_allele_data()

    elif args.rebuild:
        create_conversion_tables()
        create_fasta_ref_from_hla_allele_data()

    elif args.version:
        if args.version not in versions:
            sys.exit("[reference] Error: invalid version.")
        checkout_imgt_hla_db(versions[args.version])
        create_fasta_ref_from_hla_allele_data()
        create_conversion_tables()

    elif args.commit:
        checkout_imgt_hla_db(args.commit)
        create_fasta_ref_from_hla_allele_data()
        create_conversion_tables()

    else:
        ensure_ref_exists()
        hla_dat_commit(True)

    hline()
    log.info("")
# -------------------------------------------------------------------------------
