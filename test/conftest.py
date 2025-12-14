"""
Common fixtures, etc.
"""

import pytest

from pathlib import Path


@pytest.fixture(scope="session")
def repo_root() -> str:
    """
    Provide the path to the repository root as defined as the parent to the dir in which
    this file lives.
    """
    return str(Path(__file__).parent.parent)
