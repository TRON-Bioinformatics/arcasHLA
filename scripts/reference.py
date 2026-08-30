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

import argparse
import hashlib
import io
import logging as log
import numpy as np

import os
import shutil
import subprocess
import tarfile

from argparse import RawTextHelpFormatter
from scipy import stats
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import config
from arcas_utilities import *
from ref_paths import CORE_REFERENCE_FILES, REFERENCE_SCHEMA

# -------------------------------------------------------------------------------
#   Paths and filenames
# -------------------------------------------------------------------------------

ROOT_DIR = Path(__file__).resolve().parent.parent
REQUIRED_IMGT_FILES = ("hla.dat", "wmda/hla_nom_g.txt", "wmda/hla_nom_p.txt")
CUSTOM_FASTAS = (
    "GRCh38.all.noHLA.fasta",
    "GRCh38.chr6.noHLA.fasta",
    "GRCh38.chr6.HLA.fasta",
)


@dataclass(frozen=True)
class ReferencePaths:
    imgt_dir: Path
    out_dir: Path
    hla_dat: Path
    hla_nom_g: Path
    hla_nom_p: Path
    hla_convert_json: Path
    hla_fa: Path
    partial_fa: Path
    hla_json: Path
    partial_json: Path
    hla_idx: Path
    partial_idx: Path
    parameters_json: Path
    manifest_json: Path


def paths_for(imgt_dir, out_dir):
    """Return all source and output paths for a reference build."""
    imgt_dir = Path(imgt_dir).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    return ReferencePaths(
        imgt_dir=imgt_dir,
        out_dir=out_dir,
        hla_dat=imgt_dir / "hla.dat",
        hla_nom_g=imgt_dir / "wmda" / "hla_nom_g.txt",
        hla_nom_p=imgt_dir / "wmda" / "hla_nom_p.txt",
        hla_convert_json=out_dir / "ref" / "hla.convert.json",
        hla_fa=out_dir / "ref" / "hla.fasta",
        partial_fa=out_dir / "ref" / "hla_partial.fasta",
        hla_json=out_dir / "ref" / "hla.p.json",
        partial_json=out_dir / "ref" / "hla_partial.p.json",
        hla_idx=out_dir / "ref" / "hla.idx",
        partial_idx=out_dir / "ref" / "hla_partial.idx",
        parameters_json=out_dir / "info" / "parameters.json",
        manifest_json=out_dir / "manifest.json",
    )


# -------------------------------------------------------------------------------
#   Fetch and process IMGTHLA database
# -------------------------------------------------------------------------------


def get_mode(lengths):
    return stats.mode(lengths)[0]


def validate_imgt_source(imgt_dir):
    """Validate a read-only IMGT/HLA source tree."""
    source = Path(imgt_dir).expanduser().resolve()
    if not source.is_dir():
        sys.exit(f"[reference] Error: IMGT/HLA source is not a directory: {source}")
    missing = [
        relative
        for relative in REQUIRED_IMGT_FILES
        if not (source / relative).is_file()
    ]
    if missing:
        sys.exit(
            "[reference] Error: IMGT/HLA source is missing required files: "
            + ", ".join(missing)
        )
    return source


def _read_git_commit(imgt_dir):
    if not (Path(imgt_dir) / ".git").exists():
        return None
    result = subprocess.run(
        ["git", "-C", str(imgt_dir), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _detect_release_from_files(imgt_dir):
    candidates = (
        Path(imgt_dir) / "release_version.txt",
        Path(imgt_dir) / "Allelelist.txt",
        Path(imgt_dir) / "wmda" / "hla_nom_g.txt",
    )
    pattern = re.compile(r"(?<!\d)(\d+\.\d+\.\d+)(?!\d)")
    for candidate in candidates:
        if not candidate.is_file():
            continue
        with candidate.open("r", encoding="UTF-8", errors="replace") as file:
            for _, line in zip(range(100), file):
                match = pattern.search(line)
                if match:
                    return match.group(1)
    return None


def detect_imgt_version(imgt_dir, version_override=None):
    """Detect the IMGT release label and optional git commit read-only."""
    commit = _read_git_commit(imgt_dir)
    detected = None
    if commit:
        detected = next(
            (
                version
                for version, known_commit in config.versions.items()
                if known_commit == commit
            ),
            None,
        )
    if not detected:
        detected = _detect_release_from_files(imgt_dir)
    version = version_override or detected
    if not version:
        sys.exit(
            "[reference] Error: unable to detect the IMGT/HLA release; "
            "pass `--version X.Y.Z`."
        )
    return version, commit


def select_alleles(alleles, genes=None, max_alleles_per_gene=None):
    """Deterministically restrict alleles to a subset of genes and a maximum
    number of alleles per gene. Alleles are chosen in sorted order so that
    the same source always yields the same selection.
    """
    if not genes and not max_alleles_per_gene:
        return set(alleles)

    wanted = {gene.upper() for gene in genes} if genes else None
    selected = set()
    counts = defaultdict(int)
    for allele in sorted(alleles):
        gene = get_gene(allele)
        if wanted is not None and gene not in wanted:
            continue
        if max_alleles_per_gene and counts[gene] >= max_alleles_per_gene:
            continue
        counts[gene] += 1
        selected.add(allele)
    return selected


def process_hla_dat(hla_dat_path, genes=None, max_alleles_per_gene=None):
    """Processes IMGTHLA database, returning HLA sequences, exon locations,
    lists of complete and partial alleles and possible exon combinations.

    `genes` and `max_alleles_per_gene` build a reduced reference covering
    only part of the database.
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

    with open(hla_dat_path, "r", encoding="UTF-8") as file:
        lines = file.read().splitlines()

    # Check if hla.dat failed to download
    if len(lines) < 10:
        sys.exit(f"[reference] Error: {hla_dat_path} empty or corrupted.")

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
        elif line.startswith("FT") and re.search(r'allele\="HLA-', line):
            allele = re.split("HLA-", re.sub(r'["\n]', "", line))[1]
            gene = get_gene(allele)

            exon = sequence = False
            seq = ""

        # Exon coordinates
        elif line.startswith("FT") and re.search("exon", line):
            info = re.split(r"\s+", line)
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
        elif line.startswith("FT") and (re.search(r"\sUTR\s", line)):
            info = re.split(r"\s+", line)
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

    complete_alleles = select_alleles(complete_alleles, genes, max_alleles_per_gene)
    partial_alleles = select_alleles(partial_alleles, genes, max_alleles_per_gene)

    if genes or max_alleles_per_gene:
        if not complete_alleles and not partial_alleles:
            sys.exit(
                "[reference] Error: the requested gene selection matched no "
                "alleles; available genes are " + ", ".join(sorted(gene_set)) + "."
            )
        selected = complete_alleles | partial_alleles
        sequences = {
            allele: sequence
            for allele, sequence in sequences.items()
            if allele in selected
        }
        exons = {
            allele: coords for allele, coords in exons.items() if allele in selected
        }
        utrs = {allele: coords for allele, coords in utrs.items() if allele in selected}
        gene_set = {get_gene(allele) for allele in selected}
        gene_exons = {
            gene: numbers for gene, numbers in gene_exons.items() if gene in gene_set
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


def process_hla_nom(hla_nom):
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


def build_allele_groups(hla_nom_g, allele_keys, genes=None):
    """Build 2-field allele equivalence sets used by customize --grouping."""
    wanted = {gene.upper() for gene in genes} if genes else None
    groups = {allele: {allele} for allele in allele_keys}
    with Path(hla_nom_g).open("r", encoding="UTF-8") as file:
        for line in file:
            if line.startswith("#"):
                continue
            gene, alleles, group = line.rstrip("\n").split(";")
            if not group:
                continue
            if wanted is not None and gene.rstrip("*").upper() not in wanted:
                continue
            members = {
                process_allele(gene + allele, 2) for allele in alleles.split("/")
            }
            for allele in members:
                groups.setdefault(allele, {allele}).update(members)
    return {key: sorted(value) for key, value in sorted(groups.items())}


def _run_kallisto_index(fasta, index, reference_type, jobs=1):
    command = ["kallisto", "index", "-i", str(index)]
    if jobs and int(jobs) > 1:
        command.extend(["-t", str(jobs)])
    command.append(str(fasta))
    log.info(
        "[reference] indexing %s reference with Kallisto:\n\n\t%s\n",
        reference_type,
        " ".join(command),
    )
    return subprocess.run(command, capture_output=True, text=True, check=False)


def write_reference(
    sequences,
    info,
    fasta,
    idx,
    database,
    reference_type,
    version_token,
    jobs=1,
):
    """Writes and idxes HLA references."""
    with open(fasta, "w") as file:
        SeqIO.write(sequences, file, "fasta")

    with open(database, "w") as file:
        if len(info) == 4:
            json.dump(
                [
                    version_token,
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
                    version_token,
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

    result = _run_kallisto_index(fasta, idx, reference_type, jobs)
    if result.returncode != 0:
        raise RuntimeError(
            f"kallisto failed to index {reference_type} reference: "
            f"{result.stderr.strip()}"
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


def build_fasta(
    paths,
    version_token,
    jobs=1,
    keep_going=False,
    genes=None,
    max_alleles_per_gene=None,
):
    """Constructs HLA reference from processed sequences and exon locations."""

    log.info("[reference] IMGT/HLA database version:\n\n%s", version_token)

    log.info("[reference] Processing IMGT/HLA database")

    # Constructs cDNA sequences for alleles and adds UTRs to the set of
    # non-coding sequences
    def build_complete(allele):
        gene = get_gene(allele)
        allele_exons = sorted(exons[allele].items())
        coords = [[start, stop] for n, (start, stop) in allele_exons]
        seq = [sequences[allele][start:stop] for start, stop in coords]
        seq = "".join(seq)

        exon, exon_length = final_exon_length.get(gene, (None, None))
        if exon is not None and exon in exons[allele]:
            start, stop = exons[allele][exon]
            if stop - start > exon_length:
                seq = seq[: exon_length - (stop - start) + 1]

        cDNA[seq].add(allele)
        cdna_by_allele[process_allele(allele, 2)].add(seq)
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
    ) = process_hla_dat(paths.hla_dat, genes, max_alleles_per_gene)

    exon_combinations = get_exon_combinations()

    gene_length = defaultdict(list)
    cDNA = defaultdict(set)
    cdna_by_allele = defaultdict(set)
    combo = {str(i): defaultdict(set) for i in exon_combinations}
    other = set()

    # Build sequences for each allele
    for allele in complete_alleles:
        build_complete(allele)
        build_combination(allele)

    for allele in partial_alleles:
        build_combination(allele)

    cDNA = {seq: sorted(alleles) for seq, alleles in cDNA.items()}
    # Sorted so that UTR record order, and therefore the record identifiers
    # and file checksums, do not depend on string hash randomization.
    other = sorted(other)
    gene_length = {g: get_mode(lengths) for g, lengths in gene_length.items()}

    errors = []

    def write_or_collect(*args):
        try:
            write_reference(*args)
        except RuntimeError as error:
            if not keep_going:
                raise
            errors.append(str(error))

    log.info("[reference] Building HLA database")
    seq_out, allele_idx, lengths = complete_records(cDNA, other)
    write_or_collect(
        seq_out,
        [gene_set, allele_idx, lengths, gene_length],
        paths.hla_fa,
        paths.hla_idx,
        paths.hla_json,
        "complete",
        version_token,
        jobs,
    )

    log.info("[reference] Building partial HLA database")
    seq_out, allele_idx, lengths, exon_idx = partial_records(combo, other)
    partial_exons = {allele: exons[allele] for allele in partial_alleles}
    write_or_collect(
        seq_out,
        [gene_set, allele_idx, exon_idx, lengths, partial_exons, partial_alleles],
        paths.partial_fa,
        paths.partial_idx,
        paths.partial_json,
        "partial",
        version_token,
        jobs,
    )

    cdna = {
        allele: sorted(values, key=lambda sequence: (len(sequence), sequence))
        for allele, values in sorted(cdna_by_allele.items())
    }
    cdna_single = {allele: values[0] for allele, values in cdna.items()}
    allele_groups = build_allele_groups(paths.hla_nom_g, cdna, genes)
    return cdna, cdna_single, allele_groups, errors


def filter_convert(allele_to_group, genes=None):
    """Restricts conversion tables to a subset of genes."""
    if not genes:
        return allele_to_group

    wanted = {gene.upper() for gene in genes}
    return {
        fields: {
            allele: group
            for allele, group in mapping.items()
            if get_gene(allele) in wanted
        }
        for fields, mapping in allele_to_group.items()
    }


def build_convert(paths, genes=None):
    """Creates conversion tables for arcasHLA convert."""

    log.info("[reference] Building nomenclature conversion tables.")

    p_group = filter_convert(process_hla_nom(paths.hla_nom_p), genes)
    g_group = filter_convert(process_hla_nom(paths.hla_nom_g), genes)

    with open(paths.hla_convert_json, "w") as file:
        json.dump([p_group, g_group], file)


class NumpyEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.int64):
            return int(o)
        return json.JSONEncoder.default(self, o)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_or_link(source, destination):
    source = Path(source)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == destination.resolve():
        return
    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def _stage_custom_fastas(static_ref_dir, output_ref_dir):
    archive_path = static_ref_dir / "customization_reference_fastas.tar.gz"
    missing = []
    for name in CUSTOM_FASTAS:
        source = static_ref_dir / name
        destination = output_ref_dir / name
        if source.is_file():
            _copy_or_link(source, destination)
        else:
            missing.append(name)
    if not missing:
        return
    if not archive_path.is_file():
        raise RuntimeError(
            "missing customization FASTAs and archive: " + ", ".join(missing)
        )
    with tarfile.open(archive_path, "r:*") as archive:
        names = set(archive.getnames())
        absent = [name for name in missing if name not in names]
        if absent:
            raise RuntimeError("customization archive is missing: " + ", ".join(absent))
        for name in missing:
            source = archive.extractfile(archive.getmember(name))
            if source is None:
                raise RuntimeError(f"unable to read {name} from customization archive")
            destination = output_ref_dir / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("wb") as file:
                shutil.copyfileobj(source, file)


def _build_hla_transcripts(static_ref_dir):
    fasta_path = static_ref_dir / "GRCh38.chr6.HLA.fasta"
    if fasta_path.is_file():
        records = SeqIO.parse(fasta_path, "fasta")
    else:
        archive_path = static_ref_dir / "customization_reference_fastas.tar.gz"
        if not archive_path.is_file():
            return None
        archive = tarfile.open(archive_path, "r:*")
        source = archive.extractfile("GRCh38.chr6.HLA.fasta")
        if source is None:
            archive.close()
            return None
        records = SeqIO.parse(io.TextIOWrapper(source, encoding="UTF-8"), "fasta")

    transcripts = defaultdict(list)
    try:
        for record in records:
            match = re.search(r"gene_symbol:HLA-(\S+)", record.description)
            if match:
                transcripts[match.group(1)].append(record.id)
    finally:
        if not fasta_path.is_file():
            archive.close()
    return dict(transcripts)


def _kallisto_version():
    for command in (["kallisto", "version"], ["kallisto", "--version"]):
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode == 0:
            output = (result.stdout or result.stderr).strip()
            if output:
                return output
    return None


class ReferenceBuilder:
    """Builds a self-contained, versioned arcasHLA reference directory from a
    read-only IMGT/HLA source. Static, version-independent arcasHLA assets
    (parameters, priors, decoys, and the customization FASTAs) are staged
    from `static_data_dir` (the repository's `dat/` by default).
    """

    def __init__(
        self,
        imgt_dir,
        out_dir,
        version=None,
        force=False,
        jobs=1,
        keep_going=False,
        skip_customize=False,
        genes=None,
        max_alleles_per_gene=None,
        static_data_dir=ROOT_DIR / "dat",
    ):
        self.paths = paths_for(imgt_dir, out_dir)
        self.version_override = version
        self.force = force
        self.jobs = jobs
        self.keep_going = keep_going
        self.skip_customize = skip_customize
        self.genes = sorted({gene.upper() for gene in genes}) if genes else None
        self.max_alleles_per_gene = max_alleles_per_gene
        self.static_data_dir = Path(static_data_dir)

    def build(self):
        validate_imgt_source(self.paths.imgt_dir)
        if self.paths.manifest_json.exists() and not self.force:
            sys.exit(
                f"[reference] Error: {self.paths.out_dir} already contains "
                "manifest.json; use --force or choose another --outdir."
            )

        self.paths.out_dir.mkdir(parents=True, exist_ok=True)
        (self.paths.out_dir / "ref").mkdir(exist_ok=True)
        (self.paths.out_dir / "info").mkdir(exist_ok=True)
        if self.paths.manifest_json.exists():
            self.paths.manifest_json.unlink()

        version, commit = detect_imgt_version(
            self.paths.imgt_dir, self.version_override
        )
        version_token = commit or version

        for name in ("parameters.json", "hla_freq.tsv", "decoys_alts.json"):
            source = self.static_data_dir / "info" / name
            if not source.is_file():
                raise RuntimeError(f"missing static arcasHLA asset: {source}")
            _copy_or_link(source, self.paths.out_dir / "info" / name)

        cdna, cdna_single, allele_groups, errors = build_fasta(
            self.paths,
            version_token=version_token,
            jobs=self.jobs,
            keep_going=self.keep_going,
            genes=self.genes,
            max_alleles_per_gene=self.max_alleles_per_gene,
        )
        build_convert(self.paths, self.genes)

        tables = {
            "cDNA.json": cdna,
            "cDNA.single.json": cdna_single,
            "allele_groups.json": allele_groups,
        }
        for name, table in tables.items():
            with (self.paths.out_dir / "ref" / name).open(
                "w", encoding="UTF-8"
            ) as file:
                json.dump(table, file, sort_keys=True)

        hla_transcripts = _build_hla_transcripts(self.static_data_dir / "ref")
        if hla_transcripts is None:
            transcript_source = self.static_data_dir / "ref" / "hla_transcripts.json"
            if not transcript_source.is_file():
                raise RuntimeError(
                    f"missing static arcasHLA asset: {transcript_source}"
                )
            with transcript_source.open("r", encoding="UTF-8") as file:
                hla_transcripts = json.load(file)
        with (self.paths.out_dir / "ref" / "hla_transcripts.json").open(
            "w", encoding="UTF-8"
        ) as file:
            json.dump(hla_transcripts, file, sort_keys=True)

        if self.skip_customize:
            for name in CUSTOM_FASTAS:
                destination = self.paths.out_dir / "ref" / name
                if destination.exists():
                    destination.unlink()
        else:
            _stage_custom_fastas(
                self.static_data_dir / "ref", self.paths.out_dir / "ref"
            )

        if errors:
            raise RuntimeError("\n".join(errors))

        missing = [
            relative
            for relative in CORE_REFERENCE_FILES
            if not (self.paths.out_dir / relative).is_file()
        ]
        if missing:
            raise RuntimeError(
                "reference build did not produce required files: " + ", ".join(missing)
            )

        produced_files = set(CORE_REFERENCE_FILES)
        produced_files.update(
            {
                "ref/allele_groups.json",
                "ref/cDNA.json",
                "ref/cDNA.single.json",
                "ref/hla_transcripts.json",
            }
        )
        if not self.skip_customize:
            produced_files.update(f"ref/{name}" for name in CUSTOM_FASTAS)
        files = {
            relative: {"sha256": _sha256(self.paths.out_dir / relative)}
            for relative in sorted(produced_files)
        }

        manifest = {
            "arcashla_ref_schema": REFERENCE_SCHEMA,
            "built_at": datetime.now(timezone.utc)
            .replace(microsecond=0)
            .isoformat()
            .replace("+00:00", "Z"),
            "imgt_source": str(self.paths.imgt_dir),
            "imgt_version": version,
            "imgt_commit": commit,
            "kallisto_version": _kallisto_version(),
            "selection": {
                "genes": self.genes or "all",
                "max_alleles_per_gene": self.max_alleles_per_gene,
                "minified": bool(self.genes or self.max_alleles_per_gene),
            },
            "customization_assets": (
                "omitted"
                if self.skip_customize
                else "staged from arcasHLA static assets"
            ),
            "files": files,
        }
        with self.paths.manifest_json.open("w", encoding="UTF-8") as file:
            json.dump(manifest, file, indent=2, sort_keys=True)
            file.write("\n")
        log.info("[reference] Built reference at %s", self.paths.out_dir)
        return manifest


# -------------------------------------------------------------------------------
#   Main
# -------------------------------------------------------------------------------


def do_reference_build(
    imgt,
    outdir,
    version=None,
    force=False,
    jobs=1,
    keep_going=False,
    skip_customize=False,
    genes=None,
    max_alleles_per_gene=None,
    verbose=False,
):
    if jobs < 1:
        sys.exit("[reference] Error: --jobs must be at least 1.")
    if max_alleles_per_gene is not None and max_alleles_per_gene < 1:
        sys.exit("[reference] Error: --max-alleles-per-gene must be at least 1.")
    if verbose:
        log.basicConfig(level=log.DEBUG, format="%(message)s")
    else:
        log.basicConfig(format="%(message)s")

    try:
        return ReferenceBuilder(
            imgt,
            outdir,
            version=version,
            force=force,
            jobs=jobs,
            keep_going=keep_going,
            skip_customize=skip_customize,
            genes=genes,
            max_alleles_per_gene=max_alleles_per_gene,
        ).build()
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        sys.exit(f"[reference] Error: {error}")


def arg_check_genes(parser, argument):
    """Validates a comma separated list of HLA genes."""
    genes = {gene.strip().upper() for gene in argument.split(",") if gene.strip()}
    unknown = sorted(genes - config.genes)
    if unknown:
        parser.error(
            "invalid gene(s): "
            + ", ".join(unknown)
            + "\noptions: "
            + ", ".join(sorted(config.genes))
        )
    return sorted(genes)


def build_arg_parser(super_parser=None, subcommand_name="reference"):
    parser_args = {
        "prog": "arcasHLA reference",
        "usage": "%(prog)s build [options]",
        "add_help": False,
        "formatter_class": RawTextHelpFormatter,
    }

    if not super_parser:
        parser = argparse.ArgumentParser(**parser_args)

    else:
        parser = super_parser.add_parser(name=subcommand_name, **parser_args)

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="show this help message and exit\n\n",
        default=argparse.SUPPRESS,
    )

    ref_subparsers = parser.add_subparsers(
        dest="ref_command",
        required=True,
        metavar="build",
        help="reference subcommand\n\n",
    )

    build_parser = ref_subparsers.add_parser(
        "build",
        prog="arcasHLA reference build",
        usage="%(prog)s --imgt path --outdir path [options]",
        add_help=False,
        formatter_class=RawTextHelpFormatter,
    )

    build_parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="show this help message and exit\n\n",
        default=argparse.SUPPRESS,
    )

    build_parser.add_argument(
        "--imgt",
        required=True,
        help="read-only IMGT/HLA source directory\n\n",
        metavar="",
    )

    build_parser.add_argument(
        "--outdir",
        required=True,
        help="reference output directory\n\n",
        metavar="",
    )

    build_parser.add_argument(
        "--version",
        help="override the detected release label\n\n",
        default=None,
        metavar="",
    )

    build_parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite a previously built reference\n\n",
        default=False,
    )

    build_parser.add_argument(
        "--jobs",
        type=int,
        help="Kallisto indexing jobs\n  default: 1\n\n",
        default=1,
        metavar="",
    )

    build_parser.add_argument(
        "--keep-going",
        action="store_true",
        help="attempt both indexes before reporting failures\n\n",
        default=False,
    )

    build_parser.add_argument(
        "--skip-customize",
        action="store_true",
        help="omit large GRCh38 assets used by customize\n\n",
        default=False,
    )

    build_parser.add_argument(
        "--genes",
        help="comma separated list of HLA genes to include\n"
        + "  default: all genes in the IMGT/HLA source\n\n",
        default=None,
        type=lambda argument: arg_check_genes(build_parser, argument),
        metavar="",
    )

    build_parser.add_argument(
        "--max-alleles-per-gene",
        type=int,
        help="include at most this many alleles per gene\n\n",
        default=None,
        metavar="",
    )

    build_parser.add_argument("-v", "--verbose", action="count", default=0)

    build_parser.set_defaults(
        run_function=lambda parsed_args: do_reference_build(
            parsed_args.imgt,
            parsed_args.outdir,
            parsed_args.version,
            parsed_args.force,
            parsed_args.jobs,
            parsed_args.keep_going,
            parsed_args.skip_customize,
            parsed_args.genes,
            parsed_args.max_alleles_per_gene,
            parsed_args.verbose,
        )
    )

    return parser


def main(args):
    parser = build_arg_parser()
    parsed_args = parser.parse_args(args)

    parsed_args.run_function(parsed_args)


if __name__ == "__main__":
    main(sys.argv[1:])

# -------------------------------------------------------------------------------
