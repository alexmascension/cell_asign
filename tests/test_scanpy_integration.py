"""Optional integration tests requiring scanpy.

These tests are marked with ``scanpy`` so they can be run separately with:

    pytest -q -m scanpy

The core CI can either skip them automatically when scanpy is not installed,
or exclude them explicitly with:

    pytest -q -m "not scanpy"
"""

import pytest

from cellassign import assign_cats


@pytest.mark.scanpy
def test_assign_cats_runs_on_scanpy_pbmc3k_processed_basic():
    sc = pytest.importorskip("scanpy")

    adata = sc.datasets.pbmc3k_processed()
    adata = adata[:500].copy()

    candidate_markers = {
        "T_cells": ["CD3D", "CD3E", "TRAC", "IL7R"],
        "B_cells": ["MS4A1", "CD79A"],
        "Myeloid": ["LYZ", "S100A8", "FCGR3A", "LST1"],
    }
    markers = {
        category: [gene for gene in genes if gene in adata.var_names]
        for category, genes in candidate_markers.items()
    }
    markers = {category: genes for category, genes in markers.items() if genes}

    if len(markers) < 2:
        pytest.skip("Not enough candidate PBMC marker genes are present in this Scanpy dataset.")

    if "connectivities" not in adata.obsp:
        sc.pp.neighbors(adata)
    adata.uns.setdefault("neighbors", {}).setdefault("params", {}).setdefault("n_neighbors", 15)

    scores = assign_cats(
        adata,
        markers,
        column_groupby="louvain",
        key_added="cell_type_marker",
        min_score=0.0,
        do_return=True,
        verbose=False,
        use_raw=True,
    )

    assert "cell_type_marker" in adata.obs
    assert adata.obs["cell_type_marker"].notna().all()
    assert set(scores.columns) == set(markers)
    assert scores.shape[0] == adata.obs["louvain"].nunique()
    summary = adata.uns["cell_assign"]["cell_type_marker"]["cluster_summary"]
    assert summary["loo_stability"].between(0, 1).all()
