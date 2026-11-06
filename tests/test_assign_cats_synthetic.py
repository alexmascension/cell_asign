"""Fast synthetic tests for marker-based category assignment.

The toy matrices have no non-marker expressed genes, so the null model is
disabled (``null_correction=False``); it is tested in test_assign_cats_methodology.

These tests intentionally avoid scanpy. They construct the minimal AnnData
fields used by ``assign_cats`` directly:

- ``adata.X``
- ``adata.obs[column_groupby]``
- ``adata.var_names``
- ``adata.obsp["connectivities"]``
- ``adata.uns["neighbors"]["params"]["n_neighbors"]``
"""

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from cellassign import assign_cats

MARKERS = {
    "T_cells": ["T_marker"],
    "B_cells": ["B_marker"],
}


def make_synthetic_adata(*, as_sparse: bool = True) -> AnnData:
    """Create a tiny two-cluster AnnData object with clear marker structure."""
    x = np.array(
        [
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 5.0, 0.0],
        ]
    )
    obs = pd.DataFrame(
        {"leiden": ["0", "0", "0", "1", "1", "1"]},
        index=[f"cell_{i}" for i in range(x.shape[0])],
    )
    var = pd.DataFrame(index=["T_marker", "B_marker", "noise"])

    adata = AnnData(csr_matrix(x) if as_sparse else x, obs=obs, var=var)

    # Identity connectivities make the expected behaviour deterministic while
    # still exercising the same graph-smoothing code path.
    adata.obsp["connectivities"] = csr_matrix(np.eye(adata.n_obs))
    adata.uns["neighbors"] = {"params": {"n_neighbors": 1}}
    return adata


def test_assign_cats_assigns_expected_labels_on_sparse_input():
    adata = make_synthetic_adata(as_sparse=True)

    result = assign_cats(
        adata,
        MARKERS,
        column_groupby="leiden",
        key_added="cell_type",
        min_score=0.2,
        verbose=False,
        null_correction=False,
    )

    assert result is None
    assert adata.obs.loc[adata.obs["leiden"] == "0", "cell_type"].unique().tolist() == ["T_cells"]
    assert adata.obs.loc[adata.obs["leiden"] == "1", "cell_type"].unique().tolist() == ["B_cells"]
    assert "cell_assign" in adata.uns
    assert "cell_type" in adata.uns["cell_assign"]


def test_assign_cats_assigns_expected_labels_on_dense_input():
    adata = make_synthetic_adata(as_sparse=False)

    assign_cats(
        adata,
        MARKERS,
        column_groupby="leiden",
        key_added="cell_type",
        min_score=0.2,
        verbose=False,
        null_correction=False,
    )

    assert adata.obs.loc[adata.obs["leiden"] == "0", "cell_type"].unique().tolist() == ["T_cells"]
    assert adata.obs.loc[adata.obs["leiden"] == "1", "cell_type"].unique().tolist() == ["B_cells"]


def test_assign_cats_do_return_returns_cluster_scores():
    adata = make_synthetic_adata()

    scores = assign_cats(
        adata,
        MARKERS,
        column_groupby="leiden",
        key_added="cell_type",
        min_score=0.2,
        do_return=True,
        verbose=False,
        null_correction=False,
    )

    assert list(scores.index) == ["0", "1"]
    assert list(scores.columns) == ["T_cells", "B_cells"]
    assert scores.loc["0", "T_cells"] > scores.loc["0", "B_cells"]
    assert scores.loc["1", "B_cells"] > scores.loc["1", "T_cells"]


def test_assign_cats_respects_min_score_unassigned():
    adata = make_synthetic_adata()

    assign_cats(
        adata,
        MARKERS,
        column_groupby="leiden",
        key_added="cell_type",
        min_score=1.1,
        others_name="unassigned",
        verbose=False,
        null_correction=False,
    )

    assert set(adata.obs["cell_type"]) == {"unassigned"}


def test_assign_cats_handles_missing_marker_gene_without_failing():
    adata = make_synthetic_adata()
    markers_with_missing = {
        "T_cells": ["T_marker", "NOT_A_GENE"],
        "B_cells": ["B_marker"],
    }

    assign_cats(
        adata,
        markers_with_missing,
        column_groupby="leiden",
        key_added="cell_type",
        min_score=0.2,
        verbose=False,
        null_correction=False,
    )

    assert "cell_type" in adata.obs
    assert set(adata.obs["cell_type"]) == {"T_cells", "B_cells"}

    # The optimized implementation records missing genes in metadata; the
    # legacy implementation only prints them. Keep this assertion conditional
    # so the test suite remains backward-compatible during the refactor.
    metadata = adata.uns["cell_assign"]["cell_type"]
    if "missing_genes" in metadata:
        assert metadata["missing_genes"] == {"T_cells": ["NOT_A_GENE"]}


def test_assign_cats_intermediate_states_for_tied_categories():
    adata = make_synthetic_adata()
    tied_markers = {
        "type_a": ["T_marker"],
        "type_b": ["T_marker"],
    }

    assign_cats(
        adata,
        tied_markers,
        column_groupby="leiden",
        key_added="cell_type",
        min_score=0.2,
        intermediate_states=True,
        diff=0.01,
        verbose=False,
        null_correction=False,
    )

    assert set(adata.obs.loc[adata.obs["leiden"] == "0", "cell_type"]) == {"type_a/type_b"}
