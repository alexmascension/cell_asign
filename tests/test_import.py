from __future__ import annotations

import cellassign


def test_package_imports() -> None:
    assert cellassign.__version__ == "1.0.1"
    assert callable(cellassign.assign_cats)
