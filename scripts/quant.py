#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------
#   quant.py: genotypes from extracted chromosome 6 reads.
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

import os
import re
import json
import pickle
import sys
import argparse

import pandas as pd

from argparse import RawTextHelpFormatter
from collections import defaultdict

from arcas_utilities import *


def do_quantification(
    file,
    sample=None,
    ref=None,
    avg=200,
    std=20,
    single=False,
    LOH=False,
    purity=1.0,
    ploidy=2.0,
    threads="1",
    outdir="./",
    keep_files=False,
    temp="/tmp/",
    verbose=False,
):
    paired = not single

    if sample == None:
        sample = os.path.basename(file[0]).split(".")[0]
    else:
        sample = sample

    outdir = check_path(outdir)
    temp = create_temp(temp)

    indv_idx = ref + ".idx"
    indv_p = ref + ".p"
    indv_abundance = outdir + sample + ".quant.tsv"
    allele_results_json = outdir + sample + ".quant.alleles.json"
    gene_results_json = outdir + sample + ".quant.genes.json"
    allele_results_tsv = outdir + sample + ".quant.alleles.tsv"
    gene_results_tsv = outdir + sample + ".quant.genes.tsv"
    loh_results_tsv = outdir + sample + ".quant.loh.tsv"

    with open(indv_p, "rb") as json_file:
        genes, genotype, hla_idx, allele_idx, lengths = pickle.load(json_file)

    idx_allele = defaultdict(set)
    for idx, gene in allele_idx.items():
        idx_allele[gene].add(idx)

    if file[0].endswith(".fq.gz") or file[0].endswith(".fastq.gz"):

        command = ["kallisto", "quant", "-i", indv_idx, "-o", temp, "-t", threads]

        if single:
            command.extend(["--single -l", str(avg), "-s", str(std)])

        command.extend(file)

        output = run_command(
            command, "[quant] Quantifying with Kallisto: "
        ).stderr.decode()

        if verbose:
            print(output)

        total_reads = re.findall("(?<=processed ).+(?= reads,)", output)[0]
        total_reads = int(re.sub(",", "", total_reads))
        aligned_reads = re.findall("(?<=reads, ).+(?= reads pseudoaligned)", output)[0]
        aligned_reads = int(re.sub(",", "", aligned_reads))

        run_command(["mv", temp + "/abundance.tsv", indv_abundance])
        kallisto_results = pd.read_csv(indv_abundance, sep="\t")

    else:
        with open(file[1], "r") as json_file:
            previous_results = json.load(json_file)

        kallisto_results = pd.read_csv(file[0], sep="\t")

    idx_allele = defaultdict(set)
    hla_indices = set()
    for idx, gene in allele_idx.items():
        if gene[:-1] in genes:
            idx_allele[gene].add(idx)
            hla_indices.add(int(idx))

    lengths = defaultdict(float)
    counts = defaultdict(float)
    tpm = defaultdict(float)
    for gene, indices in idx_allele.items():
        for idx in indices:
            counts[gene] += kallisto_results.loc[int(idx)]["est_counts"]
            lengths[gene] += kallisto_results.loc[int(idx)]["length"]
            tpm[gene] += kallisto_results.loc[int(idx)]["tpm"]

    gene_results = {gene: defaultdict(int) for gene in genes}

    allele_results = {gene: defaultdict(float) for gene in genes}

    total_hla_count = 0
    for allele_id, allele in genotype.items():
        allele_results[allele_id[:-1]]["allele" + allele_id[-1]] = allele
        total_hla_count += counts[allele_id]

    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            gene_results[gene]["count"] += round(counts[allele_id])
            gene_results[gene]["tpm"] += round(tpm[allele_id])
            if counts[allele_id]:
                gene_results[gene]["abundance"] += counts[allele_id] / total_hla_count

            allele_results[gene]["allele" + allele_id[-1] + "_count"] = round(
                counts[allele_id]
            )
            allele_results[gene]["allele" + allele_id[-1] + "_tpm"] = round(
                tpm[allele_id]
            )
    for gene, allele_ids in genes.items():
        for allele_id in set(allele_ids):
            baf = allele_results[gene]["allele1_count"] / (
                allele_results[gene]["allele1_count"]
                + allele_results[gene]["allele2_count"]
            )
            allele_results[gene]["baf"] = round(min(baf, 1 - baf), 2)
    for gene in genes:
        gene_results[gene]["abundance"] = (
            str(round(gene_results[gene]["abundance"] * 100, 2)) + "%"
        )

    df = pd.DataFrame(allele_results).T
    df.index.names = ["gene"]
    try:
        df = df[
            [
                "allele1",
                "allele2",
                "allele1_count",
                "allele2_count",
                "allele1_tpm",
                "allele2_tpm",
                "baf",
            ]
        ]
    except:
        df = df[["allele1", "allele1_count", "allele1_tpm"]]
    df.to_csv(allele_results_tsv, sep="\t")

    df = pd.DataFrame(gene_results).T
    df.index.names = ["gene"]
    df = df[["count", "tpm", "abundance"]]
    df.to_csv(gene_results_tsv, sep="\t")

    with open(allele_results_json, "w") as file:
        json.dump(allele_results, file)

    with open(gene_results_json, "w") as file:
        json.dump(gene_results, file)

    if not keep_files:
        run_command(["rm", "-rf", temp])

    # LOH functionality
    if LOH:
        corrections_columns = []

        for gene in genes:
            corrections_columns.append(gene + "_CN_1")
            corrections_columns.append(gene + "_CN_2")
            corrections_columns.append(gene + "_LOSS")
            corrections_columns.append(gene + "_lost")

        corrections_df = pd.DataFrame(columns=corrections_columns)

        for gene in genes:
            baf1 = allele_results[gene]["allele1_count"] / (
                allele_results[gene]["allele1_count"]
                + allele_results[gene]["allele2_count"]
            )
            baf2 = 1 - baf1

            correction1 = (2 * baf1 * (1 + purity * (ploidy - 2) / 2) + purity - 1) / (
                purity
            )
            correction2 = (2 * baf2 * (1 + purity * (ploidy - 2) / 2) + purity - 1) / (
                purity
            )

            if correction1 < correction2:
                minor = correction1
                major = correction2
            else:
                minor = correction2
                major = correction1

            corrections_df.at[0, gene + "_CN_1"] = correction1
            corrections_df.at[0, gene + "_CN_2"] = correction2

            if (correction1 < 0.5) or (correction2 < 0.5):
                corrections_df.at[0, gene + "_LOSS"] = True

                if (correction1 < 0.5) and (correction2 < 0.5):
                    corrections_df.at[0, gene + "_lost"] = ",".join(
                        allele_results[gene][["allele1", "allele2"]].tolist()
                    )

                elif correction1 < 0.5:
                    corrections_df.at[0, gene + "_lost"] = allele_results[gene][
                        "allele1"
                    ]

                else:
                    corrections_df.at[0, gene + "_lost"] = allele_results[gene][
                        "allele2"
                    ]

            else:
                corrections_df.at[0, gene + "_LOSS"] = False
                corrections_df.at[0, gene + "_lost"] = "none"

        corrections_df.to_csv(loh_results_tsv, sep="\t", index=False)


# -------------------------------------------------------------------------------
def arg_check_files(parser, arg):
    for file in arg.split():
        if not os.path.isfile(file):
            parser.error("The file %s does not exist." % file)
        elif not (
            file.endswith("alignment.p")
            or file.endswith(".fq.gz")
            or file.endswith(".fastq.gz")
            or file.endswith(".tsv")
            or file.endswith(".json")
        ):
            parser.error("The format of %s is invalid." % file)
        return arg


def main(args):
    parser = argparse.ArgumentParser(
        prog="arcasHLA quant",
        usage="%(prog)s [options] FASTQs",
        add_help=False,
        formatter_class=RawTextHelpFormatter,
    )

    parser.add_argument(
        "file",
        help="list of fastq files",
        nargs="*",
        type=lambda x: arg_check_files(parser, x),
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="show this help message and exit\n\n",
        default=argparse.SUPPRESS,
    )

    parser.add_argument("--sample", help="sample name", type=str, default=None)

    parser.add_argument(
        "--ref",
        type=str,
        help='arcasHLA quant_ref path (e.g. "/path/to/ref/sample")\n  ',
        default=None,
        metavar="",
    )

    parser.add_argument(
        "-o", "--outdir", type=str, help="out directory\n\n", default="./", metavar=""
    )

    parser.add_argument(
        "--temp", type=str, help="temp directory\n\n", default="/tmp/", metavar=""
    )

    parser.add_argument(
        "--keep_files",
        action="count",
        help="keep intermediate files\n\n",
        default=False,
    )

    parser.add_argument(
        "--single",
        action="store_true",
        help="Include flag if single-end reads. Default is paired-end.\n\n",
        default=False,
    )

    parser.add_argument(
        "-l",
        "--avg",
        type=int,
        help="Estimated average fragment length "
        + "for single-end reads\n  default: 200\n\n",
        default=200,
    )

    parser.add_argument(
        "-s",
        "--std",
        type=int,
        help="Estimated standard deviation of fragment length "
        + "for single-end reads\n  default: 20\n\n",
        default=20,
    )

    parser.add_argument(
        "--LOH",
        action="store_true",
        help="Include flag for estimated loss of heterozygosity. "
        + "Must provide purity and ploidy estimates.\n\n",
        default=False,
    )

    parser.add_argument(
        "--purity",
        type=float,
        help="Estimated purity of sample\n  default: 1.0\n\n",
        default=1.0,
    )

    parser.add_argument(
        "--ploidy",
        type=int,
        help="Estimated ploidy of sample\n  default: 2.0\n\n",
        default=2.0,
    )

    parser.add_argument("-t", "--threads", type=str, default="1", metavar="")

    parser.add_argument("-v", "--verbose", action="count", default=False)

    parsed_args = parser.parse_args(args)

    do_quantification(
        parsed_args.file,
        parsed_args.sample,
        parsed_args.ref,
        parsed_args.avg,
        parsed_args.std,
        parsed_args.single,
        parsed_args.LOH,
        parsed_args.purity,
        parsed_args.ploidy,
        parsed_args.threads,
        parsed_args.outdir,
        parsed_args.keep_files,
        parsed_args.temp,
        parsed_args.verbose,
    )


if __name__ == "__main__":
    main(sys.argv[1:])


# -----------------------------------------------------------------------------
