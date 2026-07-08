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

It also adds per-category and summary score columns to `adata.obs`.

**We recommend runnning clustering with two resolution parameters:** the first one should be high, to
obtain a fine clustering (many clusters) and assign thus the same category to several clusters at once.
If the first resolution is too coarse you might risk assigning clusters "wrong", in the sense that splitting
a large cluster into smaller clusters may assign other unexpected clusters, or unexplored clusters that
require a more thorough inspection. 

## Expected input

`assign_cats` expects an `AnnData` object with:

* marker genes present in `adata.var_names`;
* a grouping column in `adata.obs`, such as `"leiden"` or `"louvain"`;
* a precomputed neighbour graph in `adata.obsp["connectivities"]`.

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

For each marker gene, `cellassign`:

1. extracts marker expression from the `AnnData` object;
2. smooths expression over the cell-cell neighbour graph;
3. rank-normalises positive expression values;
4. averages marker scores per category;
5. aggregates cell-level scores at the group or cluster level;
6. assigns each group to the category with the highest score.

The final label is mapped back to all cells in the group.

This makes the method less sensitive to noisy expression in individual cells, while keeping the annotation logic simple and interpretable.

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
    min_score=0.6,
    others_name="unassigned",
    verbose=True,
)
```

### Parameters

`adata`
AnnData object containing expression data, cell metadata, and a neighbour graph.

`dict_cats`
Dictionary mapping category names to marker genes.

```python
{
    "T_cells": ["CD3D", "CD3E"],
    "B_cells": ["MS4A1", "CD79A"],
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
Groups with maximum score below this value are assigned to `others_name`.

`others_name`
Label used for low-confidence assignments.

`verbose`
If `True`, prints marker genes that are missing from `adata.var_names`.

## Outputs

The function modifies `adata` in place.

Main output:

```python
adata.obs[key_added]
```

Additional score columns:

```python
adata.obs[f"{key_added}_max"]
adata.obs[f"{key_added}_mean"]
adata.obs[f"{key_added}_std"]
adata.obs[f"{key_added}_CV"]
adata.obs[f"{key_added}_{category}"]
```

Run parameters and missing marker genes are stored in:

```python
adata.uns["cell_assign"][key_added]
```

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
* **Try to go with specific markers.** Broad markers may lead to assignment of novel populations
to "broad-spectrum" populations. Sometimes a marker with lower expression but which ensures not
being expressed in other cell types can be a better marker.
* **Try to use a similar number of markers.** Populations with discordant numbers of markers
(an order of magnitude, or even 2x if marker expression is broad) may lead to a bias of assignment
to the population with greater ammount of markers.
* **Results depend on `quantile_gene_sel`.** Try several values to see how robust the assignment is.
If it varies a lot, probably a few markers are leading the lack of robustness.

## Citation

If you use `cellassign`, please cite this repository:

```bibtex
@software{cellassign,
  author = {Ascensión, Alex M.},
  title = {cellassign: marker-based cell category assignment for AnnData objects},
  url = {https://github.com/alexmascension/cell_asign},
  version = {1.0.1}
}
```

> Ascensión, Alex M. (2026). *cellassign: Lightweight marker-based assignment of cell categories in AnnData objects.*. Zenodo. https://doi.org/10.5281/zenodo.21264220

A `CITATION.cff` file is also provided for citation managers and GitHub citation metadata.

## License

This project is licensed under the MIT License.
