"""
Test the individual scripts for expected output.
"""

import pytest
import subprocess
import json
import os


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session", autouse=True)
def set_reference_version():
    """
    Fetch IMGT/HLA database version 3.24.0 before test suite
    """
    reference_cmd = f"{ROOT_DIR}/arcasHLA reference --version 3.24.0"
    subprocess.run(reference_cmd.split(), check=True)


@pytest.fixture(scope="session")
def extract_reads():
    """
    Extract reads before typing tests
    """
    output_dir = "test/output"
    extract_cmd = f"{ROOT_DIR}/arcasHLA extract test/test.bam -o {output_dir} -t 8 -v"
    subprocess.run(extract_cmd.split(), check=True)

    # Provide the individual extracted reads files.
    return [
        os.path.join(output_dir, outfile)
        for outfile in sorted(os.listdir(output_dir))
        if outfile.endswith(".fq.gz") and outfile.startswith("test")
    ]


def test_whole_allele_typing(extract_reads):
    whole_typing_cmd = (
        f"{ROOT_DIR}/arcasHLA genotype {extract_reads[0]} "
        f"{extract_reads[1]} -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o test/output -t 8 -v"
    )
    subprocess.run(whole_typing_cmd.split(), check=True)

    output_file = f"{ROOT_DIR}/test/output/test.genotype.json"
    expected_output = {
        "A": ["A*01:01:01", "A*03:01:01"],
        "B": ["B*39:01:01", "B*07:02:01"],
        "C": ["C*08:01:01", "C*01:02:01"],
        "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"],
        "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
        "DQB1": ["DQB1*02:02:01", "DQB1*06:09:01"],
        "DRB1": ["DRB1*10:01:01", "DRB1*14:02:01"],
    }
    with open(output_file, "r") as f:
        output = json.load(f)

    # Convert lists in output and expected output to sets since order does not matter for test
    for key in output:
        output[key] = set(output[key])
        expected_output[key] = set(expected_output[key])
    assert output == expected_output


def test_partial_allele_typing(extract_reads):
    partial_typing_cmd = (
        f"{ROOT_DIR}/arcasHLA partial {extract_reads[0]} "
        f"{extract_reads[1]} "
        "-g A,B,C,DPB1,DQB1,DQA1,DRB1 "
        f"-G {ROOT_DIR}/test/expected_output/test.genotype.json "
        f"-o test/output -t 8 -v"
    )
    subprocess.run(partial_typing_cmd.split(), check=True)

    output_file = f"{ROOT_DIR}/test/output/test.partial_genotype.json"
    expected_output = {
        "A": ["A*01:01:01", "A*03:01:01"],
        "B": ["B*07:02:01", "B*39:39:01"],
        "C": ["C*08:01:01", "C*01:02:01"],
        "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"],
        "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
        "DQB1": ["DQB1*06:04:01", "DQB1*02:02:01"],
        "DRB1": ["DRB1*03:02:01", "DRB1*14:02:01"],
    }
    with open(output_file, "r") as f:
        output = json.load(f)

    # Convert lists in output and expected output to sets since order does not matter for test
    for key in output:
        output[key] = set(output[key])
        expected_output[key] = set(expected_output[key])
    assert output == expected_output
