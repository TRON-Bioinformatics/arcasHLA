"""
End-to-end checks for the minified reference used in quick testing.

A minified reference covers a couple of genes and a handful of alleles,
so genotyping results are not biologically meaningful; they only have to
be reproducible.
"""

import json
import os
import subprocess

import pytest

from extract import do_extraction

MINIFIED_GENES = "A,E"
MINIFIED_ALLELES_PER_GENE = 3


@pytest.fixture(scope="session")
def minified_reference(repo_root, tmp_path_factory):
    """
    Build and select a minified reference covering one classical and one
    non-classical HLA gene.
    """
    imgt_dir = os.environ.get("ARCASHLA_TEST_IMGT_DIR") or os.path.join(
        repo_root, "dat", "IMGTHLA"
    )
    if not os.path.isfile(os.path.join(imgt_dir, "hla.dat")):
        pytest.skip("no IMGT/HLA source available")

    reference_dir = str(tmp_path_factory.mktemp("arcashla-minified-reference"))
    subprocess.run(
        [
            os.path.join(repo_root, "arcasHLA"),
            "reference",
            "build",
            "--imgt",
            imgt_dir,
            "--outdir",
            reference_dir,
            "--version",
            "3.24.0",
            "--genes",
            MINIFIED_GENES,
            "--max-alleles-per-gene",
            str(MINIFIED_ALLELES_PER_GENE),
            "--skip-partial",
            "--skip-customize",
        ],
        check=True,
    )
    return reference_dir


@pytest.fixture(scope="session")
def minified_reads(tmp_path_factory, minified_reference):
    output_dir = str(tmp_path_factory.mktemp("minified_reads"))
    do_extraction("test/test.bam", outdir=output_dir, reference=minified_reference)
    return [
        os.path.join(output_dir, outfile)
        for outfile in sorted(os.listdir(output_dir))
        if outfile.endswith(".fq.gz") and outfile.startswith("test")
    ]


def test_minified_reference_covers_selected_genes(minified_reference):
    with open(os.path.join(minified_reference, "manifest.json")) as file:
        manifest = json.load(file)

    assert manifest["selection"]["genes"] == MINIFIED_GENES.split(",")
    assert manifest["selection"]["minified"] is True
    assert manifest["selection"]["partial_reference"] == "omitted"

    with open(os.path.join(minified_reference, "ref", "cDNA.json")) as file:
        cdna = json.load(file)

    genes = {allele.split("*")[0] for allele in cdna}
    assert genes == set(MINIFIED_GENES.split(","))
    assert all(
        sum(allele.startswith(gene + "*") for allele in cdna)
        <= MINIFIED_ALLELES_PER_GENE
        for gene in genes
    )


def test_minified_genotyping_is_reproducible(
    repo_root, minified_reference, minified_reads, tmp_path
):
    def genotype(outdir):
        outdir.mkdir()
        subprocess.run(
            [
                os.path.join(repo_root, "arcasHLA"),
                "genotype",
                "--ref",
                minified_reference,
                "-o",
                str(outdir),
            ]
            + minified_reads,
            check=True,
        )
        with open(os.path.join(outdir, "test.genotype.json")) as file:
            return json.load(file)

    first = genotype(tmp_path / "first")
    second = genotype(tmp_path / "second")

    assert first == second
    assert set(first) <= set(MINIFIED_GENES.split(","))


def test_minified_reference_rejects_partial_typing(
    repo_root, minified_reference, minified_reads, tmp_path
):
    result = subprocess.run(
        [
            os.path.join(repo_root, "arcasHLA"),
            "partial",
            "--ref",
            minified_reference,
            "-o",
            str(tmp_path),
        ]
        + minified_reads,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "--skip-partial" in result.stdout + result.stderr
