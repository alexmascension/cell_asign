"""Smoke tests for the public package API."""

import cellassign
from cellassign import assign_cats


def test_public_api_imports():
    assert cellassign.__version__ == "1.0.1"
    assert callable(assign_cats)
