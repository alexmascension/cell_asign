# cellassign

[![CI](https://github.com/alexmascension/cell_asign/actions/workflows/ci.yml/badge.svg)](https://github.com/alexmascension/cell_asign/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/cellassign.svg)](https://pypi.org/project/cellassign/)
[![Python versions](https://img.shields.io/pypi/pyversions/cellassign.svg)](https://pypi.org/project/cellassign/)
[![License](https://img.shields.io/pypi/l/cellassign.svg)](https://github.com/alexmascension/cell_asign/blob/master/LICENSE)
[![Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://docs.astral.sh/ruff/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21264220.svg)](https://doi.org/10.5281/zenodo.21264220)


Lightweight marker-based assignment of cell categories in AnnData objects.

`cellassign` assigns cell-group or cluster-level labels using user-defined marker gene sets. It is designed for single-cell workflows where cells have already been clustered, and where a simple marker-based annotation layer is useful.

The package works directly on `AnnData` objects and uses the precomputed neighbour graph stored in `adata.obsp["connectivities"]`.

## Installation

Install the package from PyPI:

```bash
pip install cellassign
```

For development:

```bash
git clone https://github.com/alexmascension/cell_asign.git
cd cell_asign
pip install -e ".[dev]"
```

To install optional Scanpy support for examples and integration tests:

```bash
pip install -e ".[dev,scanpy]"
```

## Quick start

```python
import cellassign as ca

markers = {
    "T_cells": ["CD3D", "CD3E", "TRAC", "IL7R"],
    "B_cells": ["MS4A1", "CD79A"],
    "Myeloid": ["LYZ", "S100A8", "FCGR3A", "LST1"],
}

ca.assign_cats(
    adata,
    markers,
    column_groupby="leiden",
    key_added="cell_type_marker",
)
```

This adds the assigned category to:

```python
adata.obs["cell_type_marker"]
```

It also adds a confidence score, the margin to the second-best category and a leave-one-marker-out
stability column to `adata.obs`. Detailed tables are stored in `adata.obsm` and `adata.uns`.

**We recommend runnning clustering with two resolution parameters:** the first one should be high, to
obtain a fine clustering (many clusters) and assign thus the same category to several clusters at once.
If the first resolution is too coarse you might risk assigning clusters "wrong", in the sense that splitting
a large cluster into smaller clusters may assign other unexpected clusters, or unexplored clusters that
require a more thorough inspection. 

## Expected input

`assign_cats` expects an `AnnData` object with:

* marker genes present in `adata.var_names`;
* a grouping column in `adata.obs`, such as `"leiden"` or `"louvain"`;
* a precomputed neighbour graph in `adata.obsp["connectivities"]` (only needed when `smooth=True`, the default);
* log-normalised expression in `adata.X`, in a layer (`layer="..."`) or in `adata.raw` (`use_raw=True`).
  Scaled data (negative values) and raw counts trigger a warning.

For example, in a Scanpy workflow:

```python
import scanpy as sc

sc.pp.neighbors(adata)
sc.tl.leiden(adata)
```

Then run:

```python
ca.assign_cats(
    adata,
    markers,
    column_groupby="leiden",
    key_added="cell_type_marker",
)
```

## How it works

`cellassign`:

1. extracts marker expression from `adata.X`, a layer or `adata.raw` (log-normalised);
2. smooths expression over the neighbour graph with the row-normalised operator `D⁻¹(W + I)`:
   each cell gets the weighted mean of itself and its neighbours (weights from `connectivities`);
3. converts each marker to a **group-balanced percentile** in [0, 1]: every cluster carries the
   same total weight, so scores do not depend on how many cells each cluster has. Cells with zero
   expression get 0; ties share their average percentile, so results do not depend on cell order;
4. computes a cell score per category: the mean percentile of its positive markers minus the
   mean percentile of its negative markers (if any). A quantile (`quantile_gene_sel`) of that
   cell score in every cluster is the **category score**;
5. assigns each cluster to the category with the highest score (exact ties: alphabetical order);
6. accepts the assignment only if its **confidence** is at least `min_score`. The confidence
   compares the category score with the same statistic computed on `n_null` random control gene
   sets, matched to each marker by expression bin (`n_bins` equal-frequency bins of mean
   log-expression):

   `confidence = (score − null) / (ceiling − null)`

   where `ceiling` is the score a perfectly specific marker would reach in that cluster. 0 means
   "what random genes of similar expression give", 1 means "perfect marker". The null does not
   change which category is chosen, only whether it is accepted.

The final label is mapped back to all cells in the group.

Two quality measures are reported per cluster: the **margin** (best − second-best category score)
and the **leave-one-marker-out stability** (fraction of single-marker removals that keep the label).

## Main function

```python
ca.assign_cats(
    adata,
    dict_cats,
    column_groupby="leiden",
    quantile_gene_sel=0.7,
    do_return=False,
    intermediate_states=False,
    diff=0.05,
    key_added="assigned_cats",
    min_score=None,  # 0.3 with null_correction=True, 0.6 otherwise
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
)
```

### Parameters

`adata`
AnnData object containing expression data, cell metadata, and a neighbour graph.

`dict_cats`
Dictionary mapping category names to marker genes. A category is either a list of positive
markers, or a dict with positive (`"pos"`) and negative (`"neg"`) markers. Negative markers are
genes the category should not express; they are useful when two categories share their positive
markers:

```python
{
    "T_cells": ["CD3D", "CD3E"],
    "B_cells": ["MS4A1", "CD79A"],
    "CD4_CTL": {"pos": ["CD4", "GZMH", "GNLY"], "neg": ["CD8A", "CD8B"]},
    "CD8_TEM": {"pos": ["CD8A", "CD8B", "GZMH", "GZMK"], "neg": ["CCR7"]},
}
```

`column_groupby`
Column in `adata.obs` containing the groups or clusters to assign.

`quantile_gene_sel`
Quantile used to aggregate cell-level category scores within each group.

`do_return`
If `True`, returns the per-group category score table.

`intermediate_states`
If `True`, groups with similarly high category scores are labelled as intermediate states, for example `"T_cells/NK_cells"`.

`diff`
Maximum score difference from the best category for inclusion in an intermediate-state label.

`key_added`
Name of the output column added to `adata.obs`.

`min_score`
Minimum confidence of the chosen category; groups below it are assigned to `others_name`.
Defaults to 0.3 with `null_correction=True` and 0.6 otherwise. Lower values label more clusters;
higher values leave more clusters of types missing from `dict_cats` unassigned (see the benchmark).

`others_name`
Label used for low-confidence assignments.

`verbose`
If `True`, prints marker genes that are missing from `adata.var_names`.

`layer` / `use_raw`
Expression source: a key of `adata.layers`, or `adata.raw`. Mutually exclusive; default `adata.X`.

`smooth`
If `True`, smooth marker expression over the neighbour graph before ranking.

`balance_groups`
If `True`, every cluster has the same total weight in the percentiles. If `False`, every cell has
the same weight (scores then depend on cluster sizes).

`null_correction`, `n_bins`, `n_null`, `ctrl_size`, `random_state`
Null model used for the confidence: `n_null` control gene sets per category, each control gene
drawn from the same expression bin as its marker (`n_bins` equal-frequency bins of mean
log-expression, at most `ctrl_size` candidates per bin). Seeded by `random_state`.

`weight_shared_markers`
If `True`, each marker is weighted by 1 / (number of categories that list it). Off by default: in
the benchmark it lowered accuracy on fine-grained scales, where categories share lineage markers.

`leave_one_out`
If `True`, each marker of every category with at least two markers is dropped in turn and the
clusters are re-assigned. `{key_added}_loo_stability` is the fraction of perturbations that keep
the label. Only the cluster aggregation is recomputed, so the overhead is small.

`obs_details`
If `True`, also write per-category cluster scores, per-cell scores and mean/std/CV across
categories to `adata.obs`. Off by default to keep `adata.obs` compact.

`cell_score_prefix`
Prefix of the per-cell score columns written with `obs_details=True`. Defaults to
`f"{key_added}_cell_"`. Colliding output names raise an error instead of overwriting columns.

## Outputs

The function modifies `adata` in place.

Main output:

```python
adata.obs[key_added]
```

Columns always added to `adata.obs`:

```python
adata.obs[f"{key_added}_score"]  # confidence of the assigned category
adata.obs[f"{key_added}_margin"]  # best - second-best category score
adata.obs[f"{key_added}_loo_stability"]  # with leave_one_out=True (default)
```

Per-cell category scores go to `adata.obsm[f"{key_added}_cell_scores"]`. With `obs_details=True`
these columns are also added to `adata.obs`:

```python
adata.obs[f"{key_added}_mean"]  # mean category score across categories
adata.obs[f"{key_added}_std"]  # std across categories
adata.obs[f"{key_added}_CV"]  # std / mean
adata.obs[f"{key_added}_{category}"]  # cluster-level category score
adata.obs[f"{key_added}_cell_{category}"]  # per-cell category score
```

Run parameters, missing marker genes, marker weights, the package version and the cluster-level
tables (`cluster_scores`, `cluster_confidence`, `null_scores`, `null_sd`, `null_pvalues`,
`group_ceiling`, `cluster_summary` and `loo_labels`) are stored in:

```python
adata.uns["cell_assign"][key_added]
```

`cluster_summary` has one row per cluster: label, confidence, best and second-best score, margin,
the empirical p-value of the best category against its control sets, leave-one-out stability and
the markers whose removal changes the label.

If `do_return=True`, the function also returns a `pandas.DataFrame` with per-group category scores.

```python
scores = ca.assign_cats(
    adata,
    markers,
    column_groupby="leiden",
    key_added="cell_type_marker",
    do_return=True,
)
```

## Example with intermediate states

```python
ca.assign_cats(
    adata,
    markers,
    column_groupby="leiden",
    key_added="cell_type_marker",
    intermediate_states=True,
    diff=0.05,
)
```

Groups with several categories close to the maximum score are labelled with joined names, such as:

```text
T_cells/NK_cells
```

## Recommendations (based on experience)

* **Do you have unknown clusters? Try lowering `min_score`**. If they continue, or are assigned
to populations that you know are incorrect, it's likely it is a new population.
* **Check `margin` and `loo_stability`.** Low values flag clusters whose label depends on a single
marker or is close to a second category; `cluster_summary["loo_flipping_markers"]` names them.
* **Categories that differ only by markers they share** (for example CD4 vs CD8 cytotoxic T cells)
are hard for any positive-marker score; add negative markers (`{"pos": [...], "neg": [...]}`),
re-cluster the compartment and use the most specific markers available.
* **Use as negative markers only genes that are not transcribed in the population.** Protein-level
definitions do not always hold for mRNA: NK cells transcribe *CD3E*, and CD14 monocytes express some
*FCGR3A*. Because scores are percentiles, a low but detected expression still gives a large penalty.
`cluster_summary["loo_flipping_markers"]` lists negative markers whose removal changes a label
(`"<category>:-neg <gene>"`).
* **Try to go with specific markers.** Broad markers may lead to assignment of novel populations
to "broad-spectrum" populations. Sometimes a marker with lower expression but which ensures not
being expressed in other cell types can be a better marker.
* **Try to use a similar number of markers.** Populations with discordant numbers of markers
(an order of magnitude, or even 2x if marker expression is broad) may lead to a bias of assignment
to the population with greater ammount of markers.
* **Results depend on `quantile_gene_sel`.** Try several values to see how robust the assignment is.
If it varies a lot, probably a few markers are leading the lack of robustness.

## Benchmark

A separate project, [`cellassign-benchmark`](https://github.com/alexmascension/cellassign-benchmark), compares this version with 1.0.1, `scanpy.tl.score_genes`
and CellTypist on PBMC CITE-seq (Hao et al. 2021) and human cortex (Siletti et al. 2023), at a
coarse scale (all cells) and a fine scale (T cells, myeloid cells and interneurons re-clustered).
It also contains the ablation of every v2 component and the calibration of `min_score`.

## Citation

If you use `cellassign`, please cite this repository:

```bibtex
@software{cellassign,
  author = {Ascensión, Alex M.},
  title = {cellassign: marker-based cell category assignment for AnnData objects},
  url = {https://github.com/alexmascension/cell_asign},
  version = {2.0.0}
}
```

> Ascensión, Alex M. (2026). *cellassign: Lightweight marker-based assignment of cell categories in AnnData objects.*. Zenodo. https://doi.org/10.5281/zenodo.21264220

A `CITATION.cff` file is also provided for citation managers and GitHub citation metadata.

## License

This project is licensed under the MIT License.
