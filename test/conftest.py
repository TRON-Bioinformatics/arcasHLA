"""
Common fixtures, etc.
"""

import os
import pytest
import subprocess

from pathlib import Path

from customize import do_customization
from extract import do_extraction


@pytest.fixture(scope="session")
def repo_root() -> str:
    """
    Provide the path to the repository root as defined as the parent to the dir in which
    this file lives.
    """
    return str(Path(__file__).parent.parent)


@pytest.fixture(scope="session", autouse=True)
def set_reference_version(repo_root):
    """
    Fetch IMGT/HLA database version 3.24.0 before test suite
    """
    reference_cmd = f"{repo_root}/arcasHLA reference --version 3.24.0"
    subprocess.run(reference_cmd.split(), check=True)


@pytest.fixture(scope="session")
def extract_reads(tmp_path_factory):
    """
    Extract reads before typing tests
    """
    output_dir = str(tmp_path_factory.mktemp("extracted_reads"))
    do_extraction("test/test.bam", outdir=output_dir)

    # Provide the individual extracted reads files.
    return [
        os.path.join(output_dir, outfile)
        for outfile in sorted(os.listdir(output_dir))
        if outfile.endswith(".fq.gz") and outfile.startswith("test")
    ]


@pytest.fixture(scope="session")
def expected_output_dir(repo_root):
    return os.path.join(repo_root, "test/expected_output")


@pytest.fixture(scope="session")
def customize_reference(expected_output_dir, tmp_path_factory):
    """
    Provide a reference customized to the expected genotypes.
    """
    subject_name = "test"

    output_dir = str(tmp_path_factory.mktemp("custom_reference"))

    genotype_result_path = os.path.join(expected_output_dir, "test.genotype.json")

    do_customization(
        genotype=genotype_result_path, outdir=output_dir, subject=subject_name
    )

    # Provide the path to the customized reference.
    return output_dir + "/" + subject_name
