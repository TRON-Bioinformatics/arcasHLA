import json
import warnings
from pathlib import Path

import pytest

import ref_paths


def make_reference(path):
    path = Path(path)
    for relative in ref_paths.CORE_REFERENCE_FILES:
        output = path / relative
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("test", encoding="UTF-8")
    (path / "manifest.json").write_text(
        json.dumps({"arcashla_ref_schema": ref_paths.REFERENCE_SCHEMA}),
        encoding="UTF-8",
    )
    return path


@pytest.fixture(autouse=True)
def reset_reference_override(monkeypatch):
    ref_paths.configure_ref_dir(None)
    monkeypatch.delenv(ref_paths.REFERENCE_ENV_VAR, raising=False)
    yield
    ref_paths.configure_ref_dir(None)


def test_cli_reference_precedes_environment(tmp_path, monkeypatch):
    cli_ref = make_reference(tmp_path / "cli")
    env_ref = make_reference(tmp_path / "env")
    monkeypatch.setenv(ref_paths.REFERENCE_ENV_VAR, str(env_ref))

    assert ref_paths.get_ref_dir(cli_ref) == str(cli_ref.resolve())


def test_environment_precedes_legacy(tmp_path, monkeypatch):
    env_ref = make_reference(tmp_path / "env")
    legacy_ref = make_reference(tmp_path / "legacy")
    monkeypatch.setattr(ref_paths, "_legacy_ref_dir", legacy_ref.resolve())
    monkeypatch.setenv(ref_paths.REFERENCE_ENV_VAR, str(env_ref))

    assert ref_paths.get_ref_dir() == str(env_ref.resolve())


def test_legacy_fallback_warns(tmp_path, monkeypatch):
    legacy_ref = make_reference(tmp_path / "legacy")
    (legacy_ref / "manifest.json").unlink()
    monkeypatch.setattr(ref_paths, "_legacy_ref_dir", legacy_ref.resolve())

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert ref_paths.get_ref_dir() == str(legacy_ref.resolve())

    assert any(item.category is DeprecationWarning for item in caught)


def test_invalid_reference_has_actionable_error(tmp_path):
    with pytest.raises(SystemExit, match="reference build --help"):
        ref_paths.get_ref_dir(tmp_path / "missing")


def test_reference_without_partial_index_is_valid(tmp_path):
    path = make_reference(tmp_path / "minified")
    for relative in ref_paths.PARTIAL_REFERENCE_FILES:
        (path / relative).unlink()

    assert not ref_paths.is_valid_ref_dir(path)

    (path / "manifest.json").write_text(
        json.dumps(
            {
                "arcashla_ref_schema": ref_paths.REFERENCE_SCHEMA,
                "selection": {"partial_reference": "omitted"},
            }
        ),
        encoding="UTF-8",
    )

    assert ref_paths.is_valid_ref_dir(path)
