"""Lightweight marker-based cell category assignment for AnnData objects."""

from importlib.metadata import PackageNotFoundError, version

from .cellassign import assign_cats

__title__ = "cellassign"

try:
    __version__ = version("cellassign")
except PackageNotFoundError:  # pragma: no cover
    # Useful when imported from an unpacked source tree.
    __version__ = "2.1.0"

__all__ = ["assign_cats", "__version__"]
