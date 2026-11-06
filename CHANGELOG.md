# Changelog

## 2.0.0

Breaking release: scores, defaults and output columns change.

Scoring

- Smoothing uses the row-normalised operator `D⁻¹(W + I)` (weighted mean of each cell and its
  neighbours). `adata.uns["neighbors"]["params"]["n_neighbors"]` is no longer required;
  `smooth=False` disables smoothing.
- Ranks are replaced by group-balanced mid-rank percentiles over all cells (`balance_groups=True`):
  every cluster carries the same weight, zeros map to 0 and ties share their average percentile.
  Scores no longer depend on cluster sizes or on cell order.
- Null model (`null_correction=True`): control gene sets matched to each marker by expression bin
  (`n_bins`, `n_null`, `ctrl_size`, `random_state`). It defines the confidence
  `(score − null) / (ceiling − null)` used with `min_score` (default 0.3; 0.6 without null). The
  category is still chosen on the raw score.
- Negative markers: a category can be `{"pos": [...], "neg": [...]}`; its cell score is the mean
  percentile of the positive markers minus the mean percentile of the negative ones. Negative
  markers get matched control genes in the null model and are perturbed in leave-one-out.
- Optional a-priori marker specificity weights (`weight_shared_markers`, off by default).
- New `layer=` and `use_raw=` arguments; warnings for scaled data and raw counts.
- Duplicated genes inside a category are counted once.

Assignment and confidence

- Deterministic tie rule: exact ties go to the alphabetically first category (with a warning);
  intermediate-state labels are joined alphabetically; `diff` is inclusive.
- `{key_added}_margin` (best − second-best score) and leave-one-marker-out stability
  (`leave_one_out=True` by default, `{key_added}_loo_stability`).
- `cluster_summary`, `cluster_scores`, `cluster_confidence`, `null_scores`, `null_sd`,
  `null_pvalues`, `group_ceiling` and `loo_labels` stored in `adata.uns["cell_assign"][key_added]`.

Outputs

- By default only `{key_added}`, `{key_added}_score`, `{key_added}_margin` and
  `{key_added}_loo_stability` are added to `adata.obs`; per-cell scores go to `adata.obsm`.
  `obs_details=True` restores per-category columns, now prefixed (`{key_added}_cell_{category}`,
  `{key_added}_{category}`) and with collision checks.
- `{key_added}_CV` (with `obs_details=True`) is std / mean (it was mean / std).
- `{key_added}_max` is replaced by `{key_added}_score`.

## 1.0.1

- Modernized packaging with `pyproject.toml`.
- Moved the package to a canonical `src/` layout.
- Declared runtime and optional development dependencies in project metadata.
- Added Ruff configuration with `line-length = 99`.
- Added a GitHub Actions CI workflow for linting, formatting, import checks, and package builds.
