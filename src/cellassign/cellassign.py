"""Core marker-based cell category assignment functions."""

from __future__ import annotations

import warnings
from importlib.metadata import PackageNotFoundError, version

import numpy as np
import pandas as pd
from scipy import sparse

#: Suffixes of the cluster-level summary columns always written to ``adata.obs``.
DEFAULT_OBS_SUFFIXES = ("score", "margin", "loo_stability")
#: Extra cluster-level summary suffixes written when ``obs_details=True``.
DETAIL_OBS_SUFFIXES = ("mean", "std", "CV")

_DEFAULT_MIN_SCORE_NULL = 0.3
_DEFAULT_MIN_SCORE_RAW = 0.6
_CHUNK_GENES = 128


# --------------------------------------------------------------------------- #
# Input validation and data access
# --------------------------------------------------------------------------- #
def _require_column(adata, column_groupby: str) -> None:
    """Raise a clear error if the grouping column is not available."""
    if column_groupby not in adata.obs:
        raise KeyError(f"column_groupby={column_groupby!r} was not found in adata.obs")


def _get_connectivities(adata):
    """Return the neighbour graph used to smooth marker expression."""
    if "connectivities" not in adata.obsp:
        raise KeyError(
            "adata.obsp['connectivities'] was not found. Compute neighbours before "
            "calling assign_cats, for example with scanpy.pp.neighbors(adata), "
            "or use smooth=False."
        )
    return adata.obsp["connectivities"]


def _get_expression_source(adata, *, layer: str | None, use_raw: bool):
    """Return ``(var_names, full_matrix, getter, label)`` for the expression source.

    ``getter(genes)`` returns the (cells x genes) expression matrix for ``genes``.
    """
    if layer is not None and use_raw:
        raise ValueError("Use either layer=... or use_raw=True, not both.")

    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        raw = adata.raw
        return raw.var_names, raw.X, (lambda genes: raw[:, genes].X), "adata.raw.X"

    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"layer={layer!r} was not found in adata.layers")
        return (
            adata.var_names,
            adata.layers[layer],
            (lambda genes: adata[:, genes].layers[layer]),
            f"adata.layers[{layer!r}]",
        )

    return adata.var_names, adata.X, (lambda genes: adata[:, genes].X), "adata.X"


def _parse_category(category, spec) -> tuple[list, list]:
    """Return ``(positive, negative)`` marker lists for one category.

    A category is either a list of positive markers or a dict with a ``"pos"``
    list and an optional ``"neg"`` list.
    """
    if isinstance(spec, str):
        raise TypeError(f"Markers of {category!r} must be a list of genes, not a string")
    if isinstance(spec, dict):
        unknown = set(spec) - {"pos", "neg"}
        if unknown:
            raise ValueError(
                f"Category {category!r} has unknown keys {sorted(unknown)}; use 'pos' and 'neg'"
            )
        positive = list(dict.fromkeys(spec.get("pos", [])))
        negative = list(dict.fromkeys(spec.get("neg", [])))
        if not positive:
            raise ValueError(f"Category {category!r} needs at least one positive marker ('pos')")
        overlap = set(positive) & set(negative)
        if overlap:
            raise ValueError(
                f"Category {category!r} lists {sorted(overlap)} as both positive and negative"
            )
        return positive, negative
    return list(dict.fromkeys(spec)), []


def _split_marker_genes(dict_cats, var_names, *, verbose: bool):
    """Collect present/missing positive and negative markers, keeping category order.

    Duplicated genes inside one category are kept only once. Returns
    ``(positive_by_category, negative_by_category, unique_present_genes, missing)``;
    ``unique_present_genes`` contains positive and negative markers.
    """
    if not dict_cats:
        raise ValueError("dict_cats must contain at least one category")

    var_name_set = set(var_names)
    positive_by_category, negative_by_category = {}, {}
    unique_present_markers = []
    seen_markers = set()
    missing_by_category = {}

    for category, spec in dict_cats.items():
        positive, negative = _parse_category(category, spec)
        missing = {}
        for kind, genes, target in (
            ("pos", positive, positive_by_category),
            ("neg", negative, negative_by_category),
        ):
            present = [g for g in genes if g in var_name_set]
            absent = [g for g in genes if g not in var_name_set]
            target[category] = present
            for gene in present:
                if gene not in seen_markers:
                    unique_present_markers.append(gene)
                    seen_markers.add(gene)
            if absent:
                missing[kind] = absent
        if missing:
            # Plain list when only positive markers are missing (backwards compatible).
            missing_by_category[category] = missing["pos"] if list(missing) == ["pos"] else missing

    if verbose:
        for category, missing in missing_by_category.items():
            items = missing.items() if isinstance(missing, dict) else [("pos", missing)]
            for kind, genes in items:
                for gene in genes:
                    suffix = ", negative" if kind == "neg" else ""
                    print(f"Gene {gene} ({category}{suffix}) is not in var_names")

    return positive_by_category, negative_by_category, unique_present_markers, missing_by_category


def _marker_weights(markers_by_category, *, weight_shared_markers: bool) -> dict:
    """A-priori specificity weight per marker: 1 / number of categories using it."""
    counts = {}
    for genes in markers_by_category.values():
        for gene in genes:
            counts[gene] = counts.get(gene, 0) + 1
    if not weight_shared_markers:
        return dict.fromkeys(counts, 1.0)
    return {gene: 1.0 / count for gene, count in counts.items()}


def _as_csr(x):
    if sparse.issparse(x):
        return sparse.csr_matrix(x, dtype=float)
    return sparse.csr_matrix(np.asarray(x, dtype=float))


def _check_expression_values(x_markers, *, source: str) -> None:
    """Warn when marker values look scaled or like raw counts."""
    values = x_markers.data
    if values.size == 0:
        return
    if values.min() < 0:
        warnings.warn(
            f"Negative values found in {source}; the data look scaled (e.g. sc.pp.scale). "
            "Values <= 0 are treated as 'not expressed'. Use log-normalised data via "
            "layer=... or use_raw=True.",
            UserWarning,
            stacklevel=3,
        )
    elif np.all(np.mod(values, 1) == 0) and values.max() > 20:
        warnings.warn(
            f"{source} looks like raw counts (integer values up to {values.max():.0f}). "
            "Library size will confound the scores; use log-normalised data.",
            UserWarning,
            stacklevel=3,
        )


# --------------------------------------------------------------------------- #
# Groups
# --------------------------------------------------------------------------- #
def _make_grouping(groups: pd.Series) -> dict:
    """Pre-compute group codes and slices (same order as ``groupby(observed=True)``)."""
    categorical = pd.Categorical(groups)
    codes = np.asarray(categorical.codes)
    present = np.unique(codes[codes >= 0])
    order = np.argsort(codes, kind="stable")
    sorted_codes = codes[order]
    starts = np.searchsorted(sorted_codes, present, side="left")
    ends = np.searchsorted(sorted_codes, present, side="right")
    return {
        "codes": codes,
        "labels": pd.Index(categorical.categories[present], name=groups.name),
        "order": order,
        "slices": list(zip(starts, ends, strict=True)),
        "sizes": ends - starts,
    }


def _cell_weights(grouping: dict, *, balance_groups: bool) -> np.ndarray:
    """Weights summing to 1 used for the percentiles.

    With ``balance_groups=True`` every group carries the same total weight, so
    scores do not depend on how many cells each group has. Cells without a group
    get weight 0.
    """
    codes = grouping["codes"]
    weights = np.zeros(codes.shape[0], dtype=float)
    valid = codes >= 0
    if balance_groups:
        size_by_code = np.zeros(codes.max() + 1)
        size_by_code[np.unique(codes[valid])] = grouping["sizes"]
        weights[valid] = 1.0 / (len(grouping["sizes"]) * size_by_code[codes[valid]])
    else:
        weights[valid] = 1.0 / valid.sum()
    return weights


def _cluster_quantile(values: np.ndarray, grouping: dict, q: float) -> np.ndarray:
    """Per-group linear-interpolation quantile of each column of (cells x k) values.

    Rounded to 12 decimals so that floating-point near-ties become exact ties.
    """
    values = np.asarray(values)
    if values.ndim == 1:
        values = values[:, None]
    out = np.empty((len(grouping["slices"]), values.shape[1]), dtype=float)
    order = grouping["order"]
    for i, (start, end) in enumerate(grouping["slices"]):
        out[i] = np.quantile(values[order[start:end]], q, axis=0)
    return np.round(out, 12)


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #
def _row_normalized_graph(connectivities):
    """Return the row-stochastic smoothing operator ``D^-1 (W + I)``.

    Any existing diagonal in ``W`` is replaced by a self-loop of weight 1, so each
    smoothed value is a weighted mean of the cell and its neighbours, with the
    cell weighted like its strongest possible neighbour.
    """
    w = sparse.csr_matrix(connectivities, dtype=float)
    w = w - sparse.diags(w.diagonal())
    w = w + sparse.identity(w.shape[0], format="csr", dtype=float)
    row_sums = np.asarray(w.sum(axis=1)).ravel()
    return sparse.diags(1.0 / row_sums) @ w


def _smooth(x_csr, smoothing_operator) -> np.ndarray:
    """Smooth (cells x genes) expression; returns a dense array."""
    if smoothing_operator is None:
        return x_csr.toarray()
    smoothed = smoothing_operator @ x_csr
    smoothed = smoothed.toarray() if sparse.issparse(smoothed) else np.asarray(smoothed, float)
    # Values that differ only by floating-point summation order must tie, otherwise
    # percentiles would depend on the order of the cells.
    return np.round(smoothed, decimals=10)


def _weighted_percentiles(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Column-wise weighted mid-rank percentile in [0, 1]; values <= 0 map to 0.

    For a positive value x the score is ``W(v < x) + W(v == x) / 2``, where ``W``
    is the total weight of cells (weights sum to 1). Ties share the same score,
    so the result does not depend on the order of the cells.
    """
    out = np.zeros(values.shape, dtype=float)
    for col in range(values.shape[1]):
        v = values[:, col]
        positive = v > 0
        if not positive.any():
            continue
        order = np.argsort(v, kind="stable")
        sorted_v = v[order]
        new_block = np.r_[True, sorted_v[1:] != sorted_v[:-1]]
        block_id = np.cumsum(new_block) - 1
        block_weight = np.bincount(block_id, weights=weights[order])
        mid = np.cumsum(block_weight) - 0.5 * block_weight
        scores = np.empty_like(v)
        scores[order] = mid[block_id]
        out[:, col] = np.where(positive, scores, 0.0)
    return np.round(out, 12)


def _score_genes(get_expression, genes, smoothing_operator, weights, *, dtype=float):
    """Smoothed weighted percentiles for ``genes`` (processed in chunks)."""
    out = np.zeros((weights.shape[0], len(genes)), dtype=dtype)
    for start in range(0, len(genes), _CHUNK_GENES):
        chunk = list(genes[start : start + _CHUNK_GENES])
        smoothed = _smooth(_as_csr(get_expression(chunk)), smoothing_operator)
        out[:, start : start + len(chunk)] = _weighted_percentiles(smoothed, weights)
    return out


def _weighted_mean(matrix, w):
    return matrix @ w / w.sum()


def _group_ceiling(grouping: dict, weights: np.ndarray, q: float) -> np.ndarray:
    """Highest attainable group quantile: a gene expressed only, and above all
    other cells, in the group. Its cells span the top ``W_g`` of the weighted
    percentile scale, so the ``q`` quantile is ``1 - (1 - q) * W_g``.
    """
    group_weight = np.array(
        [weights[grouping["order"][start:end]].sum() for start, end in grouping["slices"]]
    )
    return 1.0 - (1.0 - q) * group_weight


def _normalised_excess(raw, null, ceiling):
    """``(raw - null) / (ceiling - null)``: 0 = null expectation, 1 = perfect marker."""
    denominator = ceiling - null
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(denominator > 1e-6, (raw - null) / denominator, 0.0)
    return np.round(out, 12)


# --------------------------------------------------------------------------- #
# Null model
# --------------------------------------------------------------------------- #
def _expression_bins(full_matrix, var_names, n_bins: int) -> pd.Series:
    """Equal-frequency bins of mean (log-normalised) expression for expressed genes."""
    means = np.asarray(full_matrix.mean(axis=0)).ravel()
    means = pd.Series(means, index=pd.Index(var_names))
    expressed = means[means > 0]
    if expressed.empty:
        raise ValueError("No expressed genes available to build the null model.")
    ranks = expressed.rank(method="first").to_numpy() - 1
    n_bins = max(1, min(int(n_bins), expressed.size))
    bins = np.floor(ranks * n_bins / expressed.size).astype(int)
    return pd.Series(bins, index=expressed.index)


def _control_pools(bins, marker_genes, ctrl_size, rng) -> dict:
    """Sample up to ``ctrl_size`` non-marker control genes from each marker bin.

    If a marker's bin has no candidate genes, the nearest non-empty bin is used.
    """
    candidates = bins.drop(index=[g for g in marker_genes if g in bins.index])
    if candidates.empty:
        raise ValueError(
            "No non-marker expressed genes are available for the null model; "
            "use null_correction=False."
        )
    by_bin = {b: list(genes) for b, genes in candidates.groupby(candidates).groups.items()}
    available = np.array(sorted(by_bin))

    pools, marker_bin = {}, {}
    for gene in marker_genes:
        b = int(bins.get(gene, 0))
        b = int(available[np.argmin(np.abs(available - b))])
        marker_bin[gene] = b
        if b not in pools:
            genes = by_bin[b]
            size = min(ctrl_size, len(genes))
            pools[b] = list(rng.choice(np.array(genes, dtype=object), size, replace=False))
    return pools, marker_bin


def _draw_controls(genes, marker_bin, pool_index, n_null, rng) -> np.ndarray:
    """(n_null x n_markers) column indices of control genes, matched by bin.

    Controls are drawn once per marker gene (not per category), so categories
    with the same markers get the same null. Within a replicate, markers that
    share a bin get distinct controls when the pool is large enough.
    """
    controls = np.empty((n_null, len(genes)), dtype=int)
    by_bin = {}
    for j, gene in enumerate(genes):
        by_bin.setdefault(marker_bin[gene], []).append(j)
    for b, positions in by_bin.items():
        pool = pool_index[b]
        replace = len(pool) < len(positions)
        for r in range(n_null):
            controls[r, positions] = rng.choice(pool, len(positions), replace=replace)
    return controls


# --------------------------------------------------------------------------- #
# Assignment
# --------------------------------------------------------------------------- #
def _assign_labels(
    scores, confidence, *, intermediate_states, diff, min_score, others_name
) -> pd.Series:
    """Assign one label per cluster with a deterministic tie rule.

    The category is chosen on ``scores`` (highest marker score) and accepted when
    ``confidence`` of that category is at least ``min_score``; otherwise the
    cluster gets ``others_name``.

    Categories are evaluated in alphabetical order: exact ties go to the
    alphabetically first category, and intermediate-state labels are joined
    alphabetically (so ``"A/B"`` and ``"B/A"`` never coexist).
    """
    order = sorted(scores.columns)
    scores = scores[order]
    confidence = confidence[order]
    best = scores.idxmax(axis=1)
    best_confidence = pd.Series(
        confidence.to_numpy()[np.arange(len(best)), scores.to_numpy().argmax(axis=1)],
        index=scores.index,
    )

    if intermediate_states:
        within = scores.ge(scores.max(axis=1) - diff, axis=0)
        labels = within.apply(lambda row: "/".join(row.index[row.to_numpy()]), axis=1)
    else:
        labels = best

    return labels.where(best_confidence >= min_score, others_name).astype(str), best_confidence


def _top_two(score_per_cluster) -> tuple[pd.Series, pd.Series]:
    """Return the best and second-best category score per cluster."""
    ordered = np.sort(score_per_cluster.to_numpy(), axis=1)[:, ::-1]
    best = pd.Series(ordered[:, 0], index=score_per_cluster.index)
    if ordered.shape[1] > 1:
        second = pd.Series(ordered[:, 1], index=score_per_cluster.index)
    else:
        second = pd.Series(np.nan, index=score_per_cluster.index)
    return best, second


def _map_cluster_values_to_obs(obs_group, values_by_cluster, *, numeric: bool = False):
    """Map cluster-level values back to cells.

    Pandas may preserve Categorical dtype when mapping from a categorical
    obs column, even if the mapped values are floats. For score columns we
    explicitly coerce back to numeric.
    """
    mapped = obs_group.map(values_by_cluster)
    if numeric:
        return pd.to_numeric(mapped, errors="coerce")
    return mapped


def _output_columns(key_added, categories, *, leave_one_out, obs_details, cell_score_prefix):
    names = [key_added, f"{key_added}_score", f"{key_added}_margin"]
    if leave_one_out:
        names.append(f"{key_added}_loo_stability")
    if obs_details:
        names += [f"{key_added}_{suffix}" for suffix in DETAIL_OBS_SUFFIXES]
        names += [f"{key_added}_{category}" for category in categories]
        names += [f"{cell_score_prefix}{category}" for category in categories]
    return names


def _check_output_columns(names, column_groupby):
    """Raise if output column names collide with each other or with the groups."""
    duplicated = sorted({name for name in names if names.count(name) > 1})
    if duplicated:
        raise ValueError(
            f"Output column names collide: {duplicated}. Rename the categories or "
            "change key_added / cell_score_prefix."
        )
    if column_groupby in names:
        raise ValueError(
            f"column_groupby={column_groupby!r} would be overwritten by the outputs. "
            "Rename the categories or change key_added / cell_score_prefix."
        )


def _package_version():
    try:
        return version("cellassign")
    except PackageNotFoundError:  # pragma: no cover
        return None


def assign_cats(
    adata,
    dict_cats,
    column_groupby="leiden",
    quantile_gene_sel=0.7,
    do_return=False,
    intermediate_states=False,
    diff=0.05,
    key_added="assigned_cats",
    min_score=None,
    others_name="unassigned",
    verbose=True,
    *,
    layer=None,
    use_raw=False,
    smooth=True,
    balance_groups=True,
    null_correction=True,
    n_bins=25,
    n_null=50,
    ctrl_size=50,
    random_state=0,
    weight_shared_markers=False,
    leave_one_out=True,
    obs_details=False,
    cell_score_prefix=None,
):
    """Assign marker-based cell categories to groups in an AnnData object.

    For every marker gene, expression is smoothed over the neighbour graph and
    converted to a group-balanced percentile. The score of a category in a group
    is a quantile (``quantile_gene_sel``) of the cell score: the mean percentile of
    its positive markers minus the mean percentile of its negative markers (if
    any). Each group gets the category with the highest score.

    The chosen category is accepted if its confidence is at least ``min_score``.
    With ``null_correction=True`` (default) the confidence compares the score with
    the same statistic computed on random control gene sets matched by expression:
    ``confidence = (score - null) / (ceiling - null)``, where ``ceiling`` is the
    score a perfectly specific marker would reach in that group (0 = what random
    genes of similar expression give, 1 = perfect marker). With
    ``null_correction=False`` the confidence is the score itself.

    Parameters
    ----------
    adata
        AnnData object. Needs ``adata.obsp['connectivities']`` when ``smooth=True``.
    dict_cats
        Dictionary mapping category names to marker genes. Each value is either a
        list of positive markers or a dict ``{"pos": [...], "neg": [...]}`` with
        optional negative markers (genes the category should not express), e.g.
        ``{"CD4 CTL": {"pos": ["CD4", "GZMH", "GNLY"], "neg": ["CD8A", "CD8B"]}}``.
    column_groupby
        Column in ``adata.obs`` containing the groups/clusters to assign.
    quantile_gene_sel
        Quantile used to aggregate cell-level category scores per group.
    do_return
        If ``True``, return the per-group category score table.
    intermediate_states
        If ``True``, groups with close category scores are labelled with joined
        names in alphabetical order, for example ``"cat_a/cat_b"``.
    diff
        Maximum difference from the best score for a category to be included in
        an intermediate-state label (inclusive).
    key_added
        Column name added to ``adata.obs`` with the final assignment.
    min_score
        Minimum confidence of the chosen category; groups below it are assigned
        to ``others_name``. Defaults to 0.3 with ``null_correction=True`` and 0.6
        otherwise.
    others_name
        Label used for low-confidence groups.
    verbose
        If ``True``, print missing marker genes.
    layer
        Key in ``adata.layers`` to read expression from. Defaults to ``adata.X``.
        Expression is expected to be log-normalised (e.g. ``log1p`` of normalised
        counts).
    use_raw
        If ``True``, read expression from ``adata.raw``. Mutually exclusive with
        ``layer``.
    smooth
        If ``True``, smooth marker expression with the row-normalised operator
        ``D^-1 (W + I)``: each cell gets the weighted mean of itself and its
        neighbours.
    balance_groups
        If ``True``, percentiles weight every cell by ``1 / (n_groups * group_size)``
        so that each group contributes equally and scores do not depend on group
        sizes. If ``False``, all cells weigh the same.
    null_correction
        If ``True``, compute the confidence against control gene sets (see above).
        The null does not change which category is chosen, only whether it is
        accepted.
    n_bins
        Number of equal-frequency bins of mean log-normalised expression used to
        match control genes to markers.
    n_null
        Number of random control gene sets per category.
    ctrl_size
        Maximum number of candidate control genes sampled per expression bin.
    random_state
        Seed for the control gene sampling.
    weight_shared_markers
        If ``True``, each marker is weighted by ``1 / number of categories that
        use it`` (a-priori specificity). Off by default: in the benchmark it
        lowered accuracy on fine-grained scales.
    leave_one_out
        If ``True``, re-assign each group after dropping each marker in turn
        (positive markers of categories with at least two of them, and every
        negative marker) and report the fraction of perturbations that keep the
        label.
    obs_details
        If ``True``, also write per-category group scores, per-cell scores and
        mean/std/CV across categories to ``adata.obs``. By default only the label,
        ``{key_added}_score``, ``{key_added}_margin`` and
        ``{key_added}_loo_stability`` are written; everything else is in
        ``adata.obsm`` / ``adata.uns``.
    cell_score_prefix
        Prefix of per-cell score columns when ``obs_details=True``. Defaults to
        ``f"{key_added}_cell_"``.

    Returns
    -------
    pandas.DataFrame or None
        Per-group category score table when ``do_return=True``; otherwise
        ``None``. The AnnData object is modified in place.

    Notes
    -----
    Outputs:

    - ``adata.obs[key_added]``: label; ``{key_added}_score``: confidence of the
      chosen category; ``{key_added}_margin``: best minus second-best category
      score; ``{key_added}_loo_stability``.
    - ``adata.obsm[f"{key_added}_cell_scores"]``: per-cell category scores.
    - ``adata.uns['cell_assign'][key_added]``: parameters, ``cluster_scores``,
      ``cluster_confidence``, ``null_scores``, ``null_sd``, ``null_pvalues``,
      ``group_ceiling``, ``cluster_summary``, ``loo_labels``, ``marker_weights``
      and ``missing_genes``.

    Ties between categories are resolved alphabetically.
    """

    _require_column(adata, column_groupby)
    if min_score is None:
        min_score = _DEFAULT_MIN_SCORE_NULL if null_correction else _DEFAULT_MIN_SCORE_RAW
    if cell_score_prefix is None:
        cell_score_prefix = f"{key_added}_cell_"
    category_names = list(dict_cats.keys())
    output_columns = _output_columns(
        key_added,
        category_names,
        leave_one_out=leave_one_out,
        obs_details=obs_details,
        cell_score_prefix=cell_score_prefix,
    )
    _check_output_columns(output_columns, column_groupby)

    var_names, full_matrix, get_expression, source = _get_expression_source(
        adata, layer=layer, use_raw=use_raw
    )
    if not pd.Index(var_names).is_unique:
        raise ValueError("var_names must be unique to assign marker genes reliably")

    markers_by_category, negative_by_category, marker_genes, missing_by_category = (
        _split_marker_genes(dict_cats, var_names, verbose=verbose)
    )
    if not any(markers_by_category.values()):
        raise ValueError("None of the marker genes in dict_cats were found in var_names")
    gene_weight = _marker_weights(markers_by_category, weight_shared_markers=weight_shared_markers)

    x_markers = _as_csr(get_expression(marker_genes))
    _check_expression_values(x_markers, source=source)

    groups = adata.obs[column_groupby]
    grouping = _make_grouping(groups)
    cluster_index = grouping["labels"]
    weights = _cell_weights(grouping, balance_groups=balance_groups)
    smoothing_operator = _row_normalized_graph(_get_connectivities(adata)) if smooth else None

    ranked = _weighted_percentiles(_smooth(x_markers, smoothing_operator), weights)
    marker_col = {gene: idx for idx, gene in enumerate(marker_genes)}

    # Null model: control genes matched by expression bin.
    if null_correction:
        rng = np.random.default_rng(random_state)
        bins = _expression_bins(full_matrix, var_names, n_bins)
        pools, marker_bin = _control_pools(bins, marker_genes, ctrl_size, rng)
        control_genes = list(dict.fromkeys(g for pool in pools.values() for g in pool))
        control_col = {gene: idx for idx, gene in enumerate(control_genes)}
        pool_index = {b: np.array([control_col[g] for g in pool]) for b, pool in pools.items()}
        all_controls = _draw_controls(marker_genes, marker_bin, pool_index, n_null, rng)
        ranked_controls = _score_genes(
            get_expression, control_genes, smoothing_operator, weights, dtype=np.float32
        )

    obs_raw, null_mean, null_sd, null_pvalue, cell_scores = {}, {}, {}, {}, {}
    loo_raw, loo_null, loo_names = {}, {}, {}
    for category in category_names:
        genes = markers_by_category[category]
        negatives = negative_by_category[category]
        if not genes:
            obs_raw[category] = np.zeros(len(cluster_index))
            null_mean[category] = np.zeros(len(cluster_index))
            null_sd[category] = np.full(len(cluster_index), np.nan)
            null_pvalue[category] = np.ones(len(cluster_index))
            cell_scores[category] = np.zeros(adata.n_obs)
            continue

        # Cell score = weighted mean of positive marker percentiles minus the mean
        # percentile of the negative markers (if any).
        w = np.array([gene_weight[g] for g in genes])
        sub = ranked[:, [marker_col[g] for g in genes]]
        neg_sub = ranked[:, [marker_col[g] for g in negatives]]
        positive_part = _weighted_mean(sub, w)
        negative_part = neg_sub.mean(axis=1) if negatives else 0.0
        cell = positive_part - negative_part
        obs_raw[category] = _cluster_quantile(cell, grouping, quantile_gene_sel)[:, 0]

        # Leave-one-marker-out: every positive marker (if >= 2) and every negative marker.
        perturbed, names = [], []
        if leave_one_out:
            if len(genes) > 1:
                for j, gene in enumerate(genes):
                    reduced = _weighted_mean(np.delete(sub, j, axis=1), np.delete(w, j))
                    perturbed.append(reduced - negative_part)
                    names.append(gene)
            for k, gene in enumerate(negatives):
                rest = np.delete(neg_sub, k, axis=1)
                perturbed.append(positive_part - (rest.mean(axis=1) if rest.shape[1] else 0.0))
                names.append(f"neg {gene}")
            if perturbed:
                loo_raw[category] = _cluster_quantile(
                    np.column_stack(perturbed), grouping, quantile_gene_sel
                )
                loo_names[category] = names

        if null_correction:
            controls = all_controls[:, [marker_col[g] for g in genes]]
            neg_controls = all_controls[:, [marker_col[g] for g in negatives]]
            weighted_sum = np.zeros((adata.n_obs, n_null))
            for j in range(len(genes)):
                weighted_sum += w[j] * ranked_controls[:, controls[:, j]]
            neg_sum = np.zeros((adata.n_obs, n_null))
            for k in range(len(negatives)):
                neg_sum += ranked_controls[:, neg_controls[:, k]]
            null_negative = neg_sum / len(negatives) if negatives else 0.0
            null_positive = weighted_sum / w.sum()
            null_cells = null_positive - null_negative
            null_per_set = _cluster_quantile(null_cells, grouping, quantile_gene_sel)
            null_mean[category] = null_per_set.mean(1)
            null_sd[category] = null_per_set.std(1, ddof=1) if n_null > 1 else np.nan
            exceed = (null_per_set >= obs_raw[category][:, None]).sum(1)
            null_pvalue[category] = (1 + exceed) / (1 + n_null)

            if category in loo_raw:
                columns = []
                if len(genes) > 1:
                    for j in range(len(genes)):
                        dropped = weighted_sum - w[j] * ranked_controls[:, controls[:, j]]
                        columns.append(dropped / (w.sum() - w[j]) - null_negative)
                for k in range(len(negatives)):
                    rest = neg_sum - ranked_controls[:, neg_controls[:, k]]
                    rest = rest / (len(negatives) - 1) if len(negatives) > 1 else 0.0
                    columns.append(null_positive - rest)
                loo_null[category] = np.column_stack(
                    [_cluster_quantile(c, grouping, quantile_gene_sel).mean(1) for c in columns]
                )
        cell_scores[category] = cell

    score_raw = pd.DataFrame(obs_raw, index=cluster_index)
    ceiling = _group_ceiling(grouping, weights, quantile_gene_sel)
    if null_correction:
        null_scores = pd.DataFrame(null_mean, index=cluster_index)
        confidence = pd.DataFrame(
            _normalised_excess(score_raw.to_numpy(), null_scores.to_numpy(), ceiling[:, None]),
            index=cluster_index,
            columns=score_raw.columns,
        )
    else:
        null_scores = None
        confidence = score_raw

    assign_kwargs = {
        "intermediate_states": intermediate_states,
        "diff": diff,
        "min_score": min_score,
        "others_name": others_name,
    }
    labels, best_confidence = _assign_labels(score_raw, confidence, **assign_kwargs)

    best, second = _top_two(score_raw)
    margin = best - second
    if not intermediate_states and (margin == 0).any():
        tied = list(margin.index[margin == 0])
        warnings.warn(
            f"Exact ties between top categories in groups {tied}; resolved alphabetically.",
            UserWarning,
            stacklevel=2,
        )
    cluster_summary = pd.DataFrame(
        {
            "label": labels,
            "score": best_confidence,
            "best": best,
            "second": second,
            "margin": margin,
        }
    )
    if null_correction:
        best_category = score_raw[sorted(score_raw.columns)].idxmax(axis=1)
        pvalues = pd.DataFrame(null_pvalue, index=cluster_index)
        cluster_summary["null_pvalue"] = [pvalues.loc[i, c] for i, c in best_category.items()]

    # Leave-one-marker-out stability.
    loo_labels = None
    if leave_one_out:
        perturbed_labels = {}
        for category, raw_scores in loo_raw.items():
            for j, gene in enumerate(loo_names[category]):
                perturbed_raw = score_raw.copy()
                perturbed_raw[category] = raw_scores[:, j]
                perturbed_confidence = confidence.copy()
                if null_correction:
                    perturbed_confidence[category] = _normalised_excess(
                        raw_scores[:, j], loo_null[category][:, j], ceiling
                    )
                else:
                    perturbed_confidence[category] = raw_scores[:, j]
                perturbed_labels[f"{category}:-{gene}"] = _assign_labels(
                    perturbed_raw, perturbed_confidence, **assign_kwargs
                )[0]
        loo_labels = pd.DataFrame(perturbed_labels, index=cluster_index)
        if loo_labels.shape[1]:
            changed = loo_labels.ne(labels, axis=0)
            stability = 1.0 - changed.mean(axis=1)
            flipping = changed.apply(lambda row: ",".join(row.index[row.to_numpy()]), axis=1)
        else:
            stability = pd.Series(np.nan, index=cluster_index)
            flipping = pd.Series("", index=cluster_index)
        cluster_summary["loo_stability"] = stability
        cluster_summary["loo_flipping_markers"] = flipping

    # Outputs in obs / obsm.
    def to_obs(values):
        return _map_cluster_values_to_obs(groups, values.to_dict(), numeric=True)

    adata.obs[key_added] = _map_cluster_values_to_obs(
        groups, labels.to_dict(), numeric=False
    ).astype(str)
    adata.obs[f"{key_added}_score"] = to_obs(best_confidence)
    adata.obs[f"{key_added}_margin"] = to_obs(margin)
    if leave_one_out:
        adata.obs[f"{key_added}_loo_stability"] = to_obs(cluster_summary["loo_stability"])

    cell_scores = pd.DataFrame(cell_scores, index=adata.obs_names)
    adata.obsm[f"{key_added}_cell_scores"] = cell_scores

    if obs_details:
        mean = score_raw.mean(axis=1)
        std = score_raw.std(axis=1)
        adata.obs[f"{key_added}_mean"] = to_obs(mean)
        adata.obs[f"{key_added}_std"] = to_obs(std)
        adata.obs[f"{key_added}_CV"] = to_obs(std / mean.replace(0, np.nan))
        for category in category_names:
            adata.obs[f"{key_added}_{category}"] = to_obs(score_raw[category])
            adata.obs[f"{cell_score_prefix}{category}"] = cell_scores[category].to_numpy()

    def str_index(frame):
        frame = frame.copy()
        frame.index = frame.index.astype(str)
        return frame

    metadata = {
        "dict_cats": {
            k: {key: list(genes) for key, genes in v.items()} if isinstance(v, dict) else list(v)
            for k, v in dict_cats.items()
        },
        "negative_markers": {k: v for k, v in negative_by_category.items() if v},
        "column_groupby": column_groupby,
        "quantile_gene_sel": quantile_gene_sel,
        "intermediate_states": intermediate_states,
        "diff": diff,
        "key_added": key_added,
        "min_score": min_score,
        "others_name": others_name,
        "missing_genes": missing_by_category,
        "layer": layer if layer is not None else "",
        "use_raw": use_raw,
        "smooth": smooth,
        "balance_groups": balance_groups,
        "null_correction": null_correction,
        "n_bins": n_bins,
        "n_null": n_null,
        "ctrl_size": ctrl_size,
        "random_state": random_state if random_state is not None else -1,
        "weight_shared_markers": weight_shared_markers,
        "marker_weights": gene_weight,
        "leave_one_out": leave_one_out,
        "obs_details": obs_details,
        "cellassign_version": _package_version() or "",
        "cluster_scores": str_index(score_raw),
        "cluster_summary": str_index(cluster_summary),
        "group_ceiling": dict(zip(cluster_index.astype(str), ceiling.tolist(), strict=True)),
    }
    if null_scores is not None:
        metadata["cluster_confidence"] = str_index(confidence)
        metadata["null_scores"] = str_index(null_scores)
        metadata["null_sd"] = str_index(pd.DataFrame(null_sd, index=cluster_index))
        metadata["null_pvalues"] = str_index(pd.DataFrame(null_pvalue, index=cluster_index))
    if loo_labels is not None and loo_labels.shape[1]:
        metadata["loo_labels"] = str_index(loo_labels)

    adata.uns.setdefault("cell_assign", {})[key_added] = metadata

    if do_return:
        return score_raw
    return None
