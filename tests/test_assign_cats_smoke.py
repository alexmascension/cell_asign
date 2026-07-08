from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from cellassign import assign_cats


def test_assign_cats_smoke_synthetic_anndata() -> None:
    """Minimal non-biological smoke test for the core AnnData path.

    This is intentionally tiny. Full synthetic and scanpy-backed tests belong
    to Phase C; this only keeps CI green and catches import/API breakage.
    """
    x = np.array(
        [
            [8.0, 0.0, 1.0, 0.0],
            [7.0, 0.0, 1.0, 0.0],
            [0.0, 8.0, 0.0, 1.0],
            [0.0, 7.0, 0.0, 1.0],
        ]
    )
    adata = AnnData(
        X=x,
        obs=pd.DataFrame({"leiden": ["0", "0", "1", "1"]}),
        var=pd.DataFrame(index=["A_marker", "B_marker", "noise_1", "noise_2"]),
    )
    adata.obsp["connectivities"] = csr_matrix(np.eye(adata.n_obs))
    adata.uns["neighbors"] = {"params": {"n_neighbors": 1}}

    assign_cats(
        adata,
        {"A": ["A_marker"], "B": ["B_marker"]},
        column_groupby="leiden",
        key_added="assigned_cats",
        min_score=0.0,
        verbose=False,
    )

    assert "assigned_cats" in adata.obs
    assert set(adata.obs["assigned_cats"]) == {"A", "B"}
    assert adata.uns["cell_assign"]["assigned_cats"]["column_groupby"] == "leiden"
