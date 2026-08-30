#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------
#   ref_paths.py: resolves the arcasHLA reference directory used by all scripts.
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

import json
import os
import sys
import warnings

from pathlib import Path

# -------------------------------------------------------------------------------

REFERENCE_SCHEMA = 1
REFERENCE_ENV_VAR = "ARCASHLA_REF_DIR"
CORE_REFERENCE_FILES = (
    "ref/hla.fasta",
    "ref/hla.idx",
    "ref/hla.p.json",
    "ref/hla_partial.fasta",
    "ref/hla_partial.idx",
    "ref/hla_partial.p.json",
    "ref/hla.convert.json",
    "info/parameters.json",
    "info/hla_freq.tsv",
    "info/decoys_alts.json",
)
PARTIAL_REFERENCE_FILES = (
    "ref/hla_partial.fasta",
    "ref/hla_partial.idx",
    "ref/hla_partial.p.json",
)

_repo_root = Path(__file__).resolve().parent.parent
_legacy_ref_dir = (_repo_root / "dat").resolve()
_configured_ref_dir = None


def configure_ref_dir(path=None):
    """Set the command-line reference override for the current process."""
    global _configured_ref_dir
    _configured_ref_dir = Path(path).expanduser().resolve() if path else None


def load_manifest(ref_dir):
    """Load a reference manifest, returning None for a legacy reference."""
    manifest_path = Path(ref_dir) / "manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        with manifest_path.open("r", encoding="UTF-8") as file:
            return json.load(file)
    except (OSError, ValueError):
        return None


def required_reference_files(manifest, required_files=CORE_REFERENCE_FILES):
    """Drop the partial reference files when the manifest omits them."""
    selection = (manifest or {}).get("selection") or {}
    if selection.get("partial_reference") == "omitted":
        return tuple(
            relative
            for relative in required_files
            if relative not in PARTIAL_REFERENCE_FILES
        )
    return tuple(required_files)


def is_valid_ref_dir(path, required_files=CORE_REFERENCE_FILES):
    """Return whether path contains a usable arcasHLA reference."""
    if not path:
        return False

    ref_dir = Path(path).expanduser().resolve()
    if not ref_dir.is_dir():
        return False

    manifest = load_manifest(ref_dir)
    required_files = required_reference_files(manifest, required_files)
    if any(not (ref_dir / relative).is_file() for relative in required_files):
        return False

    if ref_dir == _legacy_ref_dir:
        return (
            manifest is None or manifest.get("arcashla_ref_schema") == REFERENCE_SCHEMA
        )
    return bool(manifest and manifest.get("arcashla_ref_schema") == REFERENCE_SCHEMA)


def _invalid_reference_message(path=None):
    location = f" at {path}" if path else ""
    return (
        f"[reference] Error: no valid arcasHLA reference found{location}. "
        "Build one with `arcasHLA reference build --help`, then pass "
        "`--ref PATH` or set ARCASHLA_REF_DIR."
    )


def get_ref_dir(cli_ref=None):
    """
    Resolve the arcasHLA reference directory.

    Precedence is an explicit CLI path, the process-level CLI override,
    ARCASHLA_REF_DIR, then the legacy repository dat directory.
    """
    candidate = cli_ref or _configured_ref_dir or os.environ.get(REFERENCE_ENV_VAR)
    if candidate:
        candidate = Path(candidate).expanduser().resolve()
        if not is_valid_ref_dir(candidate):
            sys.exit(_invalid_reference_message(candidate))
        return str(candidate)

    if is_valid_ref_dir(_legacy_ref_dir):
        warnings.warn(
            "Using the repository dat/ reference is deprecated; build an external "
            "reference and set ARCASHLA_REF_DIR.",
            DeprecationWarning,
            stacklevel=2,
        )
        return str(_legacy_ref_dir)

    sys.exit(_invalid_reference_message())


def assert_ref_dir_valid(ref_dir=None, required_files=CORE_REFERENCE_FILES):
    """Exit with an actionable message unless the reference is usable."""
    resolved = Path(ref_dir or get_ref_dir()).expanduser().resolve()
    if not is_valid_ref_dir(resolved, required_files):
        sys.exit(_invalid_reference_message(resolved))

    manifest = load_manifest(resolved)
    if manifest and manifest.get("arcashla_ref_schema") != REFERENCE_SCHEMA:
        sys.exit(
            "[reference] Error: unsupported reference schema "
            f"{manifest.get('arcashla_ref_schema')!r}; expected {REFERENCE_SCHEMA}."
        )
    return str(resolved)


def ref_path(relative, ref_dir=None):
    """Resolve a path relative to an arcasHLA reference directory."""
    return str(Path(ref_dir or get_ref_dir()) / relative)


def imgt_path(imgt_dir, relative):
    """Resolve a path relative to a user-provided IMGT/HLA source."""
    return str(Path(imgt_dir).expanduser().resolve() / relative)


# -------------------------------------------------------------------------------
