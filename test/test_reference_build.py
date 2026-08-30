import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import reference
from ref_paths import CORE_REFERENCE_FILES

HLA_DAT = """\
ID   HLA00001
FT   allele="HLA-A*01:01:01"
FT   exon            1..3
FT                   /number="1"
FT   exon            4..6
FT                   /number="2"
FT   exon            7..9
FT                   /number="3"
FT   exon            10..12
FT                   /number="4"
FT   exon            13..15
FT                   /number="5"
FT   exon            16..18
FT                   /number="6"
FT   exon            19..21
FT                   /number="7"
FT   exon            22..24
FT                   /number="8"
SQ   Sequence 24 BP;
     atgaaacccgggtttaaacccggg 24
//
ID   HLA00002
FT   allele="HLA-B*07:02:01"
FT   exon            1..3
FT                   /number="1"
FT   exon            4..6
FT                   /number="2"
FT   exon            7..9
FT                   /number="3"
FT   exon            10..12
FT                   /number="4"
FT   exon            13..15
FT                   /number="5"
FT   exon            16..18
FT                   /number="6"
FT   exon            19..21
FT                   /number="7"
FT   exon            22..24
FT                   /number="8"
SQ   Sequence 24 BP;
     atgcccgggtttaaacccgggttt 24
//
"""


def make_imgt_source(path):
    path = Path(path)
    (path / "wmda").mkdir(parents=True)
    (path / "hla.dat").write_text(HLA_DAT, encoding="UTF-8")
    nomenclature = (
        "# version: IPD-IMGT/HLA 9.9.9\n" "A*;01:01:01;01:01P\n" "B*;07:02:01;07:02P\n"
    )
    (path / "wmda" / "hla_nom_g.txt").write_text(
        nomenclature.replace("P", "G"), encoding="UTF-8"
    )
    (path / "wmda" / "hla_nom_p.txt").write_text(nomenclature, encoding="UTF-8")
    return path


def make_static_data(path):
    path = Path(path)
    (path / "info").mkdir(parents=True)
    (path / "ref").mkdir()
    (path / "info" / "parameters.json").write_text(
        json.dumps([["A", "B"], ["prior"], {}]), encoding="UTF-8"
    )
    (path / "info" / "hla_freq.tsv").write_text(
        "allele\tprior\nA*01:01\t1\n", encoding="UTF-8"
    )
    (path / "info" / "decoys_alts.json").write_text("[]", encoding="UTF-8")
    (path / "ref" / "GRCh38.chr6.HLA.fasta").write_text(
        ">ENST00000000001.1 gene_symbol:HLA-A\nATG\n",
        encoding="UTF-8",
    )
    return path


@pytest.fixture
def fake_kallisto(monkeypatch):
    def index(fasta, index, reference_type, jobs=1):
        Path(index).write_bytes(
            b"index:" + reference_type.encode() + b":" + Path(fasta).read_bytes()
        )
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(reference, "_run_kallisto_index", index)
    monkeypatch.setattr(reference, "_kallisto_version", lambda: "kallisto 0.test")


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_build_produces_manifest_and_runtime_files(tmp_path, fake_kallisto):
    source = make_imgt_source(tmp_path / "imgt")
    static = make_static_data(tmp_path / "static")
    output = tmp_path / "reference"
    source_before = {
        path.relative_to(source): file_digest(path)
        for path in source.rglob("*")
        if path.is_file()
    }

    manifest = reference.ReferenceBuilder(
        source,
        output,
        version="9.9.9",
        skip_customize=True,
        static_data_dir=static,
    ).build()

    assert manifest["imgt_version"] == "9.9.9"
    assert manifest["arcashla_ref_schema"] == 1
    assert all((output / relative).is_file() for relative in CORE_REFERENCE_FILES)
    for name in (
        "allele_groups.json",
        "cDNA.json",
        "cDNA.single.json",
        "hla_transcripts.json",
    ):
        assert (output / "ref" / name).is_file()
    assert json.loads((output / "ref" / "hla_transcripts.json").read_text()) == {
        "A": ["ENST00000000001.1"]
    }
    for relative, metadata in manifest["files"].items():
        assert metadata["sha256"] == file_digest(output / relative)
    assert json.loads((output / "manifest.json").read_text()) == manifest
    assert source_before == {
        path.relative_to(source): file_digest(path)
        for path in source.rglob("*")
        if path.is_file()
    }


def test_build_refuses_existing_manifest(tmp_path, fake_kallisto):
    source = make_imgt_source(tmp_path / "imgt")
    static = make_static_data(tmp_path / "static")
    output = tmp_path / "reference"
    output.mkdir()
    (output / "manifest.json").write_text("{}", encoding="UTF-8")

    with pytest.raises(SystemExit, match="--force"):
        reference.ReferenceBuilder(
            source,
            output,
            version="9.9.9",
            skip_customize=True,
            static_data_dir=static,
        ).build()


def test_force_rebuild_is_deterministic(tmp_path, fake_kallisto):
    source = make_imgt_source(tmp_path / "imgt")
    static = make_static_data(tmp_path / "static")
    output = tmp_path / "reference"
    builder = reference.ReferenceBuilder(
        source,
        output,
        version="9.9.9",
        skip_customize=True,
        static_data_dir=static,
    )
    first = builder.build()["files"]

    second = reference.ReferenceBuilder(
        source,
        output,
        version="9.9.9",
        force=True,
        skip_customize=True,
        static_data_dir=static,
    ).build()["files"]

    assert first == second


@pytest.mark.parametrize("option", ["--update", "--version", "--commit", "--rebuild"])
def test_removed_legacy_modes_fail_cleanly(option):
    # Legacy top-level flags (`--update`, `--rebuild`, `--commit`, and a
    # bare `--version`) are no longer accepted: `reference` now requires a
    # `build` subcommand, wired through the same `build_arg_parser`/
    # `run_function` dispatch convention as every other script.
    arguments = [option, "3.24.0"] if option == "--version" else [option]
    with pytest.raises(SystemExit):
        reference.main(arguments)


@pytest.mark.parametrize("missing", ["directory", "hla.dat"])
def test_build_rejects_invalid_imgt_source(tmp_path, missing):
    source = tmp_path / "imgt"
    if missing == "hla.dat":
        source.mkdir()

    with pytest.raises(SystemExit, match="IMGT/HLA source"):
        reference.validate_imgt_source(source)
