"""
Common fixtures, etc.
"""

import os
import pytest
import subprocess

from pathlib import Path

from extract import do_extraction

IMGT_TEST_VERSION = "3.24.0"
# Upstream IMGT/HLA tags the 3.24.0 release as ``v3.24.0-alpha``; the bare
# version string is not a valid git ref.
IMGT_TEST_TAG = "v3.24.0-alpha"


@pytest.fixture(scope="session")
def repo_root() -> str:
    """
    Provide the path to the repository root as defined as the parent to the dir in which
    this file lives.
    """
    return str(Path(__file__).parent.parent)


@pytest.fixture(scope="session")
def set_reference_version(repo_root, tmp_path_factory):
    """
    Build and select an external IMGT/HLA 3.24.0 reference.
    """
    imgt_dir = os.environ.get("ARCASHLA_TEST_IMGT_DIR")
    if not imgt_dir:
        local_imgt = Path(repo_root) / "dat" / "IMGTHLA"
        if (local_imgt / "hla.dat").is_file():
            imgt_dir = str(local_imgt)
        else:
            imgt_dir = str(tmp_path_factory.mktemp("imgt-hla"))
            subprocess.run(
                [
                    "git",
                    "clone",
                    "--depth",
                    "1",
                    "--branch",
                    IMGT_TEST_TAG,
                    "https://github.com/ANHIG/IMGTHLA.git",
                    imgt_dir,
                ],
                check=True,
            )

    reference_dir = str(tmp_path_factory.mktemp("arcashla-reference"))
    subprocess.run(
        [
            f"{repo_root}/arcasHLA",
            "reference",
            "build",
            "--imgt",
            imgt_dir,
            "--outdir",
            reference_dir,
            "--version",
            IMGT_TEST_VERSION,
        ],
        check=True,
    )
    previous = os.environ.get("ARCASHLA_REF_DIR")
    os.environ["ARCASHLA_REF_DIR"] = reference_dir
    yield reference_dir
    if previous is None:
        os.environ.pop("ARCASHLA_REF_DIR", None)
    else:
        os.environ["ARCASHLA_REF_DIR"] = previous


@pytest.fixture(scope="session")
def extract_reads(tmp_path_factory, set_reference_version):
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
