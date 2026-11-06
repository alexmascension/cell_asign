"""Tests for smoothing, percentiles, null model, weights, tie rules and outputs."""

import warnings

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy.sparse import csr_matrix

from cellassign import assign_cats
from cellassign.cellassign import (
    _cluster_quantile,
    _make_grouping,
    _row_normalized_graph,
    _weighted_percentiles,
)


def make_adata(n_a=30, n_b=30, n_background=60, seed=0):
    """Two groups; A_1/A_2 mark group 0, B_1 marks group 1, 'ubiq' is everywhere.

    Background genes are ubiquitous with a range of means, for the null model.
    """
    rng = np.random.default_rng(seed)
    n = n_a + n_b
    group = np.array(["0"] * n_a + ["1"] * n_b)
    is_a = group == "0"
    specific = np.column_stack(
        [
            np.where(is_a, rng.poisson(4, n), rng.poisson(0.05, n)),  # A_1
            np.where(is_a, rng.poisson(4, n), rng.poisson(0.05, n)),  # A_2
            np.where(~is_a, rng.poisson(4, n), rng.poisson(0.05, n)),  # B_1
            rng.poisson(2, n),  # ubiq
        ]
    )
    background = rng.poisson(rng.lognormal(-0.5, 1, n_background)[None, :], (n, n_background))
    x = np.log1p(np.column_stack([specific, background]).astype(float))
    genes = ["A_1", "A_2", "B_1", "ubiq"] + [f"bg{i}" for i in range(n_background)]
    adata = AnnData(
        csr_matrix(x),
        obs=pd.DataFrame({"leiden": group}, index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=genes),
    )
    rows, cols = [], []
    for g in ["0", "1"]:
        idx = np.where(group == g)[0]
        for i in idx:
            nbrs = rng.choice(idx[idx != i], 5, replace=False)
            rows += [i] * 5
            cols += list(nbrs)
    w = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    adata.obsp["connectivities"] = ((w + w.T) > 0).astype(float).tocsr()
    return adata


MARKERS = {"A": ["A_1", "A_2"], "B": ["B_1"]}


def labels_by_group(adata, key="ct"):
    return adata.obs.groupby("leiden", observed=True)[key].first().to_dict()


# --------------------------------------------------------------------------- #
# Building blocks
# --------------------------------------------------------------------------- #
def test_smoothing_operator_is_row_stochastic_with_self_loop():
    w = csr_matrix(np.array([[0, 1, 0.5], [1, 0, 0], [0.5, 0, 0]]))
    p = _row_normalized_graph(w).toarray()
    np.testing.assert_allclose(p.sum(axis=1), 1.0)
    np.testing.assert_allclose(np.diag(p), [1 / 2.5, 1 / 2, 1 / 1.5])
    x = np.array([4.0, 2.0, 6.0])
    assert p[0] @ x == pytest.approx((4 + 2 + 3) / 2.5)


def test_weighted_percentiles_midrank_and_zeros():
    values = np.array([[0.0], [2.0], [2.0], [5.0]])
    out = _weighted_percentiles(values, np.full(4, 0.25))[:, 0]
    np.testing.assert_allclose(out, [0.0, 0.5, 0.5, 0.875])


def test_cluster_quantile_matches_pandas():
    rng = np.random.default_rng(0)
    groups = pd.Series(pd.Categorical(rng.choice(["b", "a", "c"], 200)), name="g")
    values = rng.random((200, 3))
    expected = pd.DataFrame(values).groupby(groups, observed=True).quantile(0.7)
    got = _cluster_quantile(values, _make_grouping(groups), 0.7)
    np.testing.assert_allclose(got, expected.to_numpy())


# --------------------------------------------------------------------------- #
# Methodology
# --------------------------------------------------------------------------- #
def test_null_correction_rejects_ubiquitous_marker():
    """A non-specific gene passes the raw threshold but not the null-corrected one."""
    raw = make_adata()
    assign_cats(raw, {"Junk": ["ubiq"]}, key_added="ct", null_correction=False, verbose=False)
    assert set(labels_by_group(raw).values()) == {"Junk"}

    corrected = make_adata()
    assign_cats(corrected, {"Junk": ["ubiq"]}, key_added="ct", verbose=False)
    assert set(labels_by_group(corrected).values()) == {"unassigned"}

    good = make_adata()
    assign_cats(good, MARKERS, key_added="ct", verbose=False)
    assert labels_by_group(good) == {"0": "A", "1": "B"}


def test_normalised_score_is_one_for_perfect_marker_and_zero_for_noise():
    adata = make_adata()
    assign_cats(
        adata, {"A": ["A_1"], "Junk": ["ubiq"]}, key_added="ct", smooth=False, verbose=False
    )
    scores = adata.uns["cell_assign"]["ct"]["cluster_confidence"]
    assert scores.loc["0", "A"] > 0.8
    # With only two groups the null-normalised noise is large; it must stay below
    # the default threshold, sign aside.
    assert (scores["Junk"] < 0.3).all()


def test_null_correction_only_affects_acceptance():
    with_null, without_null = make_adata(), make_adata()
    markers = {"A": ["A_1", "A_2"], "B": ["B_1"], "Junk": ["ubiq"]}
    assign_cats(with_null, markers, key_added="ct", min_score=-np.inf, verbose=False)
    assign_cats(
        without_null,
        markers,
        key_added="ct",
        min_score=-np.inf,
        null_correction=False,
        verbose=False,
    )
    assert labels_by_group(with_null) == labels_by_group(without_null)
    summary = with_null.uns["cell_assign"]["ct"]["cluster_summary"]
    assert summary["null_pvalue"].between(0, 1).all()


def test_balanced_percentiles_do_not_depend_on_group_sizes():
    scores = {}
    for balance in (True, False):
        own = []
        for n_a in (20, 180):
            adata = make_adata(n_a=n_a, n_b=200 - n_a)
            s = assign_cats(
                adata,
                MARKERS,
                key_added="ct",
                null_correction=False,
                balance_groups=balance,
                do_return=True,
                verbose=False,
            )
            own.append(s.loc["0", "A"])
        scores[balance] = abs(own[0] - own[1])
    assert scores[True] < 0.05
    assert scores[False] > 0.2


def test_shared_markers_are_downweighted():
    adata = make_adata()
    markers = {"A": ["A_1", "ubiq"], "B": ["B_1", "ubiq"]}
    assign_cats(adata, markers, key_added="ct", weight_shared_markers=True, verbose=False)
    weights = adata.uns["cell_assign"]["ct"]["marker_weights"]
    assert weights == {"A_1": 1.0, "ubiq": 0.5, "B_1": 1.0}

    assign_cats(adata, markers, key_added="ct2", verbose=False)  # default: no weighting
    assert set(adata.uns["cell_assign"]["ct2"]["marker_weights"].values()) == {1.0}


def test_scores_are_invariant_to_cell_order():
    adata = make_adata()
    perm = np.random.default_rng(1).permutation(adata.n_obs)
    shuffled = adata[perm].copy()
    for ad in (adata, shuffled):
        assign_cats(ad, MARKERS, key_added="ct", min_score=0, verbose=False)
    pd.testing.assert_frame_equal(
        adata.obsm["ct_cell_scores"],
        shuffled.obsm["ct_cell_scores"].loc[adata.obs_names],
        rtol=1e-9,
    )
    pd.testing.assert_frame_equal(
        adata.uns["cell_assign"]["ct"]["cluster_scores"],
        shuffled.uns["cell_assign"]["ct"]["cluster_scores"],
    )


def test_layer_and_use_raw_are_used():
    adata = make_adata()
    adata.layers["lognorm"] = adata.X.copy()
    adata.raw = adata.copy()
    rng = np.random.default_rng(3)
    adata.X = csr_matrix(adata.X.toarray()[rng.permutation(adata.n_obs)])  # wrong data in X

    assign_cats(adata, MARKERS, key_added="from_layer", layer="lognorm", verbose=False)
    assign_cats(adata, MARKERS, key_added="from_raw", use_raw=True, verbose=False)
    for key in ("from_layer", "from_raw"):
        assert labels_by_group(adata, key) == {"0": "A", "1": "B"}

    with pytest.raises(ValueError, match="either layer"):
        assign_cats(adata, MARKERS, layer="lognorm", use_raw=True, verbose=False)


def test_negative_values_warn():
    adata = make_adata()
    adata.X = csr_matrix(adata.X.toarray() - 1.0)
    with pytest.warns(UserWarning, match="scaled"):
        assign_cats(adata, MARKERS, null_correction=False, min_score=0, verbose=False)


# --------------------------------------------------------------------------- #
# Outputs
# --------------------------------------------------------------------------- #
def test_default_obs_columns_are_compact():
    adata = make_adata()
    before = set(adata.obs.columns)
    assign_cats(adata, MARKERS, key_added="ct", verbose=False)
    added = set(adata.obs.columns) - before
    assert added == {"ct", "ct_score", "ct_margin", "ct_loo_stability"}
    assert list(adata.obsm["ct_cell_scores"].columns) == ["A", "B"]
    meta = adata.uns["cell_assign"]["ct"]
    expected = {"cluster_scores", "cluster_confidence", "null_scores", "null_pvalues"}
    assert expected | {"cluster_summary"} <= set(meta)


def test_obs_details_and_collisions():
    adata = make_adata()
    assign_cats(adata, MARKERS, key_added="ct", obs_details=True, verbose=False)
    expected = {"ct_mean", "ct_std", "ct_CV", "ct_A", "ct_B", "ct_cell_A", "ct_cell_B"}
    assert expected <= set(adata.obs.columns)

    with pytest.raises(ValueError, match="collide"):
        assign_cats(
            adata,
            {"score": ["A_1"], "B": ["B_1"]},
            key_added="ct",
            obs_details=True,
            verbose=False,
        )
    with pytest.raises(ValueError, match="column_groupby"):
        assign_cats(
            adata,
            {"leiden": ["A_1"], "B": ["B_1"]},
            obs_details=True,
            cell_score_prefix="",
            verbose=False,
        )


def test_cv_margin_definitions():
    adata = make_adata()
    scores = assign_cats(
        adata, MARKERS, key_added="ct", obs_details=True, do_return=True, verbose=False
    )
    row = adata.obs.loc[adata.obs.leiden == "0"].iloc[0]
    raw = adata.uns["cell_assign"]["ct"]["cluster_scores"].loc["0"]
    assert row["ct_CV"] == pytest.approx(raw.std() / raw.mean())
    s = scores.loc["0"]
    assert row["ct_margin"] == pytest.approx(s.max() - s.min())
    assert s.max() <= 1.0 + 1e-9
    confidence = adata.uns["cell_assign"]["ct"]["cluster_confidence"].loc["0"]
    assert row["ct_score"] == pytest.approx(confidence[s.idxmax()])


def test_tie_rule_does_not_depend_on_dict_order():
    labels = []
    for order in (("zeta", "alpha"), ("alpha", "zeta")):
        adata = make_adata()
        markers = {name: ["A_1"] for name in order}
        with pytest.warns(UserWarning, match="resolved alphabetically"):
            assign_cats(adata, markers, key_added="ct", min_score=-1, verbose=False)
        labels.append(adata.obs["ct"].tolist())
    assert labels[0] == labels[1]
    assert set(labels[0]) == {"alpha"}


def test_intermediate_states_alphabetical_and_inclusive_diff():
    adata = make_adata()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assign_cats(
            adata,
            {"zeta": ["A_1"], "alpha": ["A_1"], "B": ["B_1"]},
            key_added="ct",
            intermediate_states=True,
            diff=0.0,
            verbose=False,
        )
    assert labels_by_group(adata)["0"] == "alpha/zeta"


def test_leave_one_out_matches_rerun_without_marker():
    adata = make_adata()
    # Dropping A_1 makes A identical to B (exact tie, resolved alphabetically to A),
    # so group "1" flips from B to A in one of the two perturbations.
    markers = {"A": ["A_1", "B_1"], "B": ["B_1"]}
    kwargs = {
        "key_added": "ct",
        "min_score": 0,
        "null_correction": False,
        "weight_shared_markers": False,
        "verbose": False,
    }
    assign_cats(adata, markers, **kwargs)
    meta = adata.uns["cell_assign"]["ct"]
    summary, loo = meta["cluster_summary"], meta["loo_labels"]

    assert list(loo.columns) == ["A:-A_1", "A:-B_1"]  # B has one marker: not perturbed
    assert summary["label"].to_dict() == {"0": "A", "1": "B"}
    assert summary["loo_stability"].to_dict() == {"0": 1.0, "1": 0.5}
    assert summary.loc["1", "loo_flipping_markers"] == "A:-A_1"

    for column in loo.columns:
        category, gene = column.split(":-")
        reduced = {
            k: [g for g in v if not (k == category and g == gene)] for k, v in markers.items()
        }
        rerun = make_adata()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assign_cats(rerun, reduced, **kwargs)
        assert loo[column].to_dict() == labels_by_group(rerun)


def test_leave_one_out_with_null_runs_and_is_bounded():
    adata = make_adata()
    assign_cats(adata, MARKERS, key_added="ct", verbose=False)
    meta = adata.uns["cell_assign"]["ct"]
    summary, loo = meta["cluster_summary"], meta["loo_labels"]
    assert list(loo.columns) == ["A:-A_1", "A:-A_2"]
    # Both A markers are informative, so dropping either keeps the labels.
    assert summary["loo_stability"].to_dict() == {"0": 1.0, "1": 1.0}
    assert adata.obs["ct_loo_stability"].between(0, 1).all()


def test_h5ad_roundtrip(tmp_path):
    adata = make_adata()
    assign_cats(adata, MARKERS, key_added="ct", obs_details=True, verbose=False)
    path = tmp_path / "out.h5ad"
    adata.write_h5ad(path)
    from anndata import read_h5ad

    back = read_h5ad(path)
    assert "cluster_summary" in back.uns["cell_assign"]["ct"]
    assert "ct_cell_scores" in back.obsm


# --------------------------------------------------------------------------- #
# Negative markers
# --------------------------------------------------------------------------- #
# "alpha" and "beta" share their only positive marker; they differ by what they
# must NOT express: group 0 expresses A_1, group 1 expresses B_1.
NEGATIVE_MARKERS = {
    "alpha": {"pos": ["ubiq"], "neg": ["A_1"]},
    "beta": {"pos": ["ubiq"], "neg": ["B_1"]},
}
NO_NULL = {"key_added": "ct", "min_score": -np.inf, "null_correction": False, "verbose": False}


def test_negative_markers_separate_categories_with_shared_positives():
    without = make_adata()
    with pytest.warns(UserWarning, match="resolved alphabetically"):
        assign_cats(without, {"alpha": ["ubiq"], "beta": ["ubiq"]}, **NO_NULL)
    assert set(labels_by_group(without).values()) == {"alpha"}

    with_negatives = make_adata()
    assign_cats(with_negatives, NEGATIVE_MARKERS, **NO_NULL)
    assert labels_by_group(with_negatives) == {"0": "beta", "1": "alpha"}

    scores = with_negatives.uns["cell_assign"]["ct"]["cluster_scores"]
    assert scores.loc["0", "alpha"] < 0  # penalised by A_1
    assert with_negatives.uns["cell_assign"]["ct"]["negative_markers"] == {
        "alpha": ["A_1"],
        "beta": ["B_1"],
    }


def test_negative_markers_mix_with_plain_lists_and_null():
    adata = make_adata()
    markers = {"A": ["A_1", "A_2"], "B": {"pos": ["B_1"], "neg": ["A_1"]}}
    assign_cats(adata, markers, key_added="ct", verbose=False)
    assert labels_by_group(adata) == {"0": "A", "1": "B"}
    confidence = adata.uns["cell_assign"]["ct"]["cluster_confidence"]
    assert np.isfinite(confidence.to_numpy()).all()


def test_negative_marker_validation():
    adata = make_adata()
    bad = [
        ({"A": "A_1"}, TypeError),
        ({"A": {"pos": ["A_1"], "not": ["B_1"]}}, ValueError),
        ({"A": {"neg": ["B_1"]}}, ValueError),
        ({"A": {"pos": ["A_1"], "neg": ["A_1"]}}, ValueError),
    ]
    for markers, error in bad:
        with pytest.raises(error):
            assign_cats(adata, markers, verbose=False)


def test_missing_negative_markers_are_reported():
    adata = make_adata()
    markers = {"A": {"pos": ["A_1", "NOPE"], "neg": ["ALSO_NOPE"]}, "B": ["B_1", "NOPE"]}
    assign_cats(adata, markers, key_added="ct", verbose=False)
    missing = adata.uns["cell_assign"]["ct"]["missing_genes"]
    assert missing == {"A": {"pos": ["NOPE"], "neg": ["ALSO_NOPE"]}, "B": ["NOPE"]}


def test_leave_one_out_drops_negative_markers():
    adata = make_adata()
    assign_cats(adata, NEGATIVE_MARKERS, **NO_NULL)
    meta = adata.uns["cell_assign"]["ct"]
    loo, summary = meta["loo_labels"], meta["cluster_summary"]
    assert list(loo.columns) == ["alpha:-neg A_1", "beta:-neg B_1"]
    # Without its negative marker, alpha (= ubiq) beats beta (= ubiq - B_1) in group 0.
    assert loo.loc["0", "alpha:-neg A_1"] == "alpha"
    assert summary.loc["0", "loo_stability"] == pytest.approx(0.5)
    assert summary.loc["0", "loo_flipping_markers"] == "alpha:-neg A_1"

    rerun = make_adata()
    reduced = {"alpha": ["ubiq"], "beta": NEGATIVE_MARKERS["beta"]}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        assign_cats(rerun, reduced, **NO_NULL)
    assert loo["alpha:-neg A_1"].to_dict() == labels_by_group(rerun)


def test_leave_one_out_with_negatives_and_null_runs():
    adata = make_adata()
    markers = {"A": {"pos": ["A_1", "A_2"], "neg": ["B_1"]}, "B": {"pos": ["B_1"], "neg": ["A_1"]}}
    assign_cats(adata, markers, key_added="ct", verbose=False)
    loo = adata.uns["cell_assign"]["ct"]["loo_labels"]
    assert list(loo.columns) == ["A:-A_1", "A:-A_2", "A:-neg B_1", "B:-neg A_1"]
    assert adata.obs["ct_loo_stability"].between(0, 1).all()


def test_h5ad_roundtrip_with_negative_markers(tmp_path):
    adata = make_adata()
    assign_cats(adata, {"A": ["A_1"], "B": {"pos": ["B_1"], "neg": ["A_1"]}}, verbose=False)
    path = tmp_path / "neg.h5ad"
    adata.write_h5ad(path)
    from anndata import read_h5ad

    back = read_h5ad(path)
    assert back.uns["cell_assign"]["assigned_cats"]["dict_cats"]["B"]["neg"].tolist() == ["A_1"]
