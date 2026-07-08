from __future__ import annotations

import pytest

pytestmark = pytest.mark.scanpy


def test_scanpy_optional_import() -> None:
    """Optional scanpy smoke test.

    The CI can run this only in jobs where scanpy is installed. More meaningful
    scanpy.datasets.pbmc3k_processed() tests are for Phase C.
    """
    scanpy = pytest.importorskip("scanpy")
    assert scanpy.__version__
