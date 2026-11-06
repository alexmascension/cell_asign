"""Smoke tests for the public package API."""

import cellassign
from cellassign import assign_cats


def test_public_api_imports():
    assert cellassign.__version__ == "2.1.0"
    assert callable(assign_cats)
