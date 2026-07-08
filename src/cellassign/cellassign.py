"""Core marker-based cell category assignment functions."""

from __future__ import annotations

import numpy as np
from scipy import sparse


def _require_column(adata, column_groupby: str) -> None:
    """Raise a clear error if the grouping column is not available."""
    if column_groupby not in adata.obs:
        raise KeyError(f"column_groupby={column_groupby!r} was not found in adata.obs")


def _get_connectivities(adata):
    """Return the neighbour graph used to smooth marker expression."""
    if "connectivities" not in adata.obsp:
        raise KeyError(
            "adata.obsp['connectivities'] was not found. Compute neighbours before "
            "calling assign_cats, for example with scanpy.pp.neighbors(adata)."
        )
    return adata.obsp["connectivities"]


def _get_n_neighbors(adata, connectivities) -> float:
    """Return the n_neighbors as a normalization factor for the joint 
    expression used to smooth out cells with extreme expression."""
    try:
        n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
    except KeyError:
        # Older or manually-created AnnData objects may contain connectivities
        # without the scanpy metadata. Falling back to the average row sum keeps
        # the function usable for synthetic AnnData objects and preserves the
        # idea of neighbour-count normalisation.
        row_sums = np.asarray(connectivities.sum(axis=1)).ravel()
        positive_row_sums = row_sums[row_sums > 0]
        if positive_row_sums.size == 0:
            raise KeyError(
                "Could not infer n_neighbors: adata.uns['neighbors']['params']"
                "['n_neighbors'] is missing and the connectivity graph has no "
                "positive row sums."
            ) from None
        n_neighbors = float(np.median(positive_row_sums))

    if n_neighbors <= 0:
        raise ValueError("n_neighbors must be a positive value")
    return float(n_neighbors)


def _split_marker_genes(dict_cats, var_names, *, verbose: bool):
    """Collect present/missing marker genes while preserving category order."""
    if not dict_cats:
        raise ValueError("dict_cats must contain at least one category")

    var_name_set = set(var_names)
    markers_by_category = {}
    unique_present_markers = []
    seen_markers = set()
    missing_by_category = {}

    for category, genes in dict_cats.items():
        present_genes = []
        missing_genes = []

        for gene in genes:
            if gene in var_name_set:
                present_genes.append(gene)
                if gene not in seen_markers:
                    unique_present_markers.append(gene)
                    seen_markers.add(gene)
            else:
                missing_genes.append(gene)

        markers_by_category[category] = present_genes
        if missing_genes:
            missing_by_category[category] = missing_genes

    if verbose:
        for missing_genes in missing_by_category.values():
            for gene in missing_genes:
                print(f"Gene {gene} is not on the list")

    return markers_by_category, unique_present_markers, missing_by_category


def _marker_expression_matrix(adata, marker_genes):
    """Return marker expression as a 2D matrix suitable for sparse products."""
    x_markers = adata[:, marker_genes].X
    if sparse.issparse(x_markers):
        return x_markers.tocsr()
    return sparse.csr_matrix(np.asarray(x_markers))


def _smooth_marker_expression(adata, marker_genes, connectivities, 
                              n_neighbors: float) -> np.ndarray:
    """Smooth all marker genes over the neighbour graph in one matrix product."""
    if not marker_genes:
        return np.zeros((adata.n_obs, 0), dtype=float)

    x_markers = _marker_expression_matrix(adata, marker_genes)
    smoothed = connectivities @ x_markers
    smoothed = smoothed / n_neighbors

    if sparse.issparse(smoothed):
        return smoothed.toarray().astype(float, copy=False)
    return np.asarray(smoothed, dtype=float)


def _rank_positive_columns(values: np.ndarray) -> np.ndarray:
    """Rank-normalise positive values column-wise to the [0, 1] interval.

    This keeps the original intent of ``np.argsort(np.argsort(...))`` but avoids
    divide-by-zero artefacts when a marker is expressed in zero or one cell.
    """
    ranked = np.zeros(values.shape, dtype=float)

    for col_idx in range(values.shape[1]):
        positive_mask = values[:, col_idx] > 0
        if not positive_mask.any():
            continue

        positive_values = values[positive_mask, col_idx]
        positive_ranks = np.argsort(np.argsort(positive_values)).astype(float)

        max_rank = positive_ranks.max()
        if max_rank > 0:
            positive_ranks /= max_rank
        else:
            positive_ranks[:] = 1.0

        ranked[positive_mask, col_idx] = positive_ranks

    return ranked


def _category_score_matrix(
    ranked_markers: np.ndarray,
    marker_genes: list,
    markers_by_category: dict,
    n_obs: int,
) -> dict:
    """Average ranked marker scores per category."""
    marker_to_col = {gene: idx for idx, gene in enumerate(marker_genes)}
    category_scores = {}

    for category, category_markers in markers_by_category.items():
        if not category_markers:
            category_scores[category] = np.zeros(n_obs, dtype=float)
            continue

        marker_cols = [marker_to_col[gene] for gene in category_markers]
        category_scores[category] = ranked_markers[:, marker_cols].mean(axis=1)

    return category_scores


def _map_cluster_stat(adata, column_groupby: str, values_by_cluster):
    """Map a per-cluster Series back to cells."""
    return adata.obs[column_groupby].map(values_by_cluster)


def assign_cats(
    adata,
    dict_cats,
    column_groupby="leiden",
    quantile_gene_sel=0.7,
    do_return=False,
    intermediate_states=False,
    diff=0.05,
    key_added="assigned_cats",
    min_score=0.6,
    others_name="unassigned",
    verbose=True,
):
    """Assign marker-based cell categories to groups in an AnnData object.

    This function uses marker gene sets assigned to different categories so
    that clusters, for example Leiden clusters, can be assigned to one category.
    It is designed for cases where assigning each cell individually would be
    noisier than assigning a precomputed cell group or cluster.

    Parameters
    ----------
    adata
        AnnData object with a precomputed neighbour graph in
        ``adata.obsp['connectivities']``.
    dict_cats
        Dictionary mapping category names to marker genes.
    column_groupby
        Column in ``adata.obs`` containing the groups/clusters to assign.
    quantile_gene_sel
        Quantile used to aggregate cell-level category scores per group.
    do_return
        If ``True``, return the per-group score table.
    intermediate_states
        If ``True``, groups with close category scores are labelled with joined
        names, for example ``"cat_a/cat_b"``.
    diff
        Maximum difference from the best score for a category to be included in
        an intermediate-state label.
    key_added
        Column name added to ``adata.obs`` with the final assignment.
    min_score
        Groups whose maximum category score is below this value are assigned to
        ``others_name``.
    others_name
        Label used for low-confidence groups.
    verbose
        If ``True``, print missing marker genes.

    Returns
    -------
    pandas.DataFrame or None
        Per-group category score table when ``do_return=True``; otherwise
        ``None``. The AnnData object is modified in place.
    """
    _require_column(adata, column_groupby)

    if not adata.var_names.is_unique:
        raise ValueError("adata.var_names must be unique to assign marker genes reliably")

    connectivities = _get_connectivities(adata)
    n_neighbors = _get_n_neighbors(adata, connectivities)

    markers_by_category, marker_genes, missing_by_category = _split_marker_genes(
        dict_cats,
        adata.var_names,
        verbose=verbose,
    )

    smoothed_markers = _smooth_marker_expression(adata, marker_genes, connectivities, n_neighbors)
    ranked_markers = _rank_positive_columns(smoothed_markers)
    category_scores = _category_score_matrix(
        ranked_markers,
        marker_genes,
        markers_by_category,
        adata.n_obs,
    )

    category_names = list(dict_cats.keys())
    for category in category_names:
        adata.obs[category] = category_scores[category]

    score_per_cluster = (
        adata.obs[[column_groupby] + category_names]
        .groupby(column_groupby)
        .quantile(quantile_gene_sel)
    )

    cluster_std = score_per_cluster.std(axis=1)
    cluster_mean = score_per_cluster.mean(axis=1)
    cluster_max = score_per_cluster.max(axis=1)

    adata.obs[f"{key_added}_std"] = _map_cluster_stat(adata, column_groupby, cluster_std)
    adata.obs[f"{key_added}_mean"] = _map_cluster_stat(adata, column_groupby, cluster_mean)
    adata.obs[f"{key_added}_max"] = _map_cluster_stat(adata, column_groupby, cluster_max)
    adata.obs[f"{key_added}_CV"] = adata.obs[f"{key_added}_mean"] / adata.obs[
        f"{key_added}_std"
    ]

    for category in score_per_cluster.columns:
        adata.obs[f"{key_added}_{category}"] = _map_cluster_stat(
            adata,
            column_groupby,
            score_per_cluster[category],
        )

    if intermediate_states:
        names_per_cluster = []
        for cluster in score_per_cluster.index:
            scores_cluster = score_per_cluster.loc[cluster]
            scores_cluster = scores_cluster[scores_cluster > scores_cluster.max() - diff]
            names_per_cluster.append("/".join(scores_cluster.index))
        final_cat_by_cluster = dict(zip(score_per_cluster.index, names_per_cluster, strict=True))
    else:
        final_cat_by_cluster = score_per_cluster.idxmax(axis=1).to_dict()

    adata.obs[key_added] = adata.obs[column_groupby].map(final_cat_by_cluster).astype(str)
    adata.obs.loc[adata.obs[f"{key_added}_max"] < min_score, key_added] = others_name

    if "cell_assign" not in adata.uns:
        adata.uns["cell_assign"] = {}

    adata.uns["cell_assign"][key_added] = {
        "dict_cats": dict_cats,
        "column_groupby": column_groupby,
        "quantile_gene_sel": quantile_gene_sel,
        "do_return": do_return,
        "intermediate_states": intermediate_states,
        "diff": diff,
        "key_added": key_added,
        "min_score": min_score,
        "others_name": others_name,
        "missing_genes": missing_by_category,
    }

    if do_return:
        return score_per_cluster
    return None

