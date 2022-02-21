import numpy as np
from scipy.sparse import csr_matrix


def assign_cats(adata, dict_cats, column_groupby='leiden', quantile_gene_sel=0.7, do_return=False, intermediate_states=False, diff=0.05,
                key_added='assigned_cats', min_score=0.6, others_name='unassigned', verbose=True):
    """
    This functions uses a set of genes assigned to different categories so that leiden clusters can be assigned to one of these categories.
    For example, to categorize fibroblasts from pericytes, endothelial cells, or cells with high mitochondrial content.
    It could be done with each cell individually, but it is better to use clusters to discern the different categories because
    the method, although efficient, can sometimes be noisy due to the noisiness of the sc datasets.
    """

    for cat in list(dict_cats.keys()):
        mat_cat = np.zeros((len(adata), len(dict_cats[cat])), dtype=float)

        for gene_idx, gene in enumerate(dict_cats[cat]):
            try:
                if type(adata.X) == np.ndarray:
                    mat_cat[:, gene_idx] = np.asarray(np.dot(adata.obsp['connectivities'], csr_matrix(
                        adata[:, gene].X)).todense()).ravel() / adata.uns['neighbors']['params']['n_neighbors']
                else:
                    mat_cat[:, gene_idx] = np.asarray(np.dot(adata.obsp['connectivities'], adata[:, gene].X).todense(
                    )).ravel() / adata.uns['neighbors']['params']['n_neighbors']
                mat_cat[mat_cat[:, gene_idx] > 0, gene_idx] = np.argsort(
                    np.argsort(mat_cat[mat_cat[:, gene_idx] > 0, gene_idx]))
                mat_cat[:, gene_idx] /= np.max(mat_cat[:, gene_idx])
            except:
                if verbose:
                    print(f'Gene {gene} is not on the list')

        sum_mat_cat = np.asarray(mat_cat.mean(1)).ravel()
        adata.obs[cat] = sum_mat_cat

    score_per_cluster = adata.obs[[column_groupby] + list(
        dict_cats.keys())].groupby(column_groupby).quantile(quantile_gene_sel)
    max_cat_dict_std = dict(zip(score_per_cluster.std(
        1).index, score_per_cluster.std(1).values))
    adata.obs[f'{key_added}_std'] = [max_cat_dict_std[i]
                                     for i in adata.obs[column_groupby]]
    max_cat_dict_mean = dict(zip(score_per_cluster.mean(
        1).index, score_per_cluster.mean(1).values))
    adata.obs[f'{key_added}_mean'] = [max_cat_dict_mean[i]
                                      for i in adata.obs[column_groupby]]
    max_cat_dict_max = dict(zip(score_per_cluster.max(
        1).index, score_per_cluster.max(1).values))
    adata.obs[f'{key_added}_max'] = [max_cat_dict_max[i]
                                     for i in adata.obs[column_groupby]]
    adata.obs[f'{key_added}_CV'] = adata.obs[f'{key_added}_mean'] / \
        adata.obs[f'{key_added}_std']

    for cat in score_per_cluster.columns:
        max_cat_dict = dict(zip(score_per_cluster.index,
                            score_per_cluster[cat].values))
        adata.obs[f'{key_added}_{cat}'] = [max_cat_dict[i]
                                           for i in adata.obs[column_groupby]]

    # For each cluster we will identify which categories are close to the highest one, and merge their names.
    if intermediate_states:
        list_names_cats_per_cluster = []
        for cluster in score_per_cluster.index:
            scores_cluster = score_per_cluster.loc[cluster]
            scores_cluster = scores_cluster[scores_cluster >
                                            scores_cluster.max() - diff]
            list_names_cats_per_cluster.append('/'.join(scores_cluster.index))

        final_cat_dict = dict(zip(score_per_cluster.idxmax(
            axis=1).index, list_names_cats_per_cluster))
    else:
        final_cat_dict = dict(zip(score_per_cluster.idxmax(
            axis=1).index, score_per_cluster.idxmax(axis=1).values))

    adata.obs[f'{key_added}'] = [str(final_cat_dict[i])
                                 for i in adata.obs[column_groupby]]

    adata.obs[f'{key_added}'][adata.obs[f'{key_added}_max']
                              < min_score] = others_name

    if do_return:
        return score_per_cluster
