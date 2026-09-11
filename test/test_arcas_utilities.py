"""
Unit tests for the error-handling behavior of arcas_utilities.run_command.

Verifies that a non-zero exit status raises a RuntimeError whose message
carries the actual exit code and the captured stderr, that
ignore_exit_status suppresses the raise, and that the RuntimeError propagates
through wrappers such as check_path.
"""

import pytest

import arcas_utilities


def test_run_command_raises_on_nonzero_exit_ls():
    """ls on a missing path (exit status 2) must raise a RuntimeError whose
    message contains both the exit status and the captured stderr."""
    with pytest.raises(RuntimeError) as excinfo:
        arcas_utilities.run_command(["ls", "/nonexistent_dir_xyz"])

    message = str(excinfo.value)
    assert "Command exited with non-zero status 2" in message
    assert "No such file or directory" in message


def test_run_command_raises_on_nonzero_exit_cat():
    """cat on a missing file (exit status 1) must raise, reporting the
    actual exit code in the message (distinct from ls' status 2)."""
    with pytest.raises(RuntimeError, match="Command exited with non-zero status 1"):
        arcas_utilities.run_command(["cat", "/nonexistent_file_xyz"])


def test_run_command_ignore_exit_status_suppresses_raise():
    """With ignore_exit_status=True a failing command must not raise, and
    the CompletedProcess (with its non-zero returncode) is returned."""
    output = arcas_utilities.run_command(
        ["ls", "/nonexistent_dir_xyz"], ignore_exit_status=True
    )

    assert output.returncode == 2


def test_check_path_propagates_run_command_error(tmp_path):
    """check_path creates its directory via run_command; a failing mkdir
    must surface as the RuntimeError rather than being swallowed."""
    with pytest.raises(RuntimeError, match="non-zero status"):
        arcas_utilities.check_path(str(tmp_path / "missing_parent" / "child"))
