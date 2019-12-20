from adata import Adata
import scscope as DeepImpute
import pandas as pd
import phenograph
from natsort import natsorted
import pickle
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import scanpy.api as sc


def _set_default_colors_for_categorical_obs(adata, value_to_plot):
    """
    Sets the adata.uns[value_to_plot + '_colors'] using default color palettes

    Parameters
    ----------
    adata : annData object
    value_to_plot : name of a valid categorical observation

    Returns
    -------
    None
    """
    from scanpy.plotting import palettes
    from matplotlib import rcParams
    from scanpy import logging as logg


    categories = adata.obs[value_to_plot].cat.categories
    length = len(categories)

    # check if default matplotlib palette has enough colors
    if len(rcParams['axes.prop_cycle'].by_key()['color']) >= length:
        cc = rcParams['axes.prop_cycle']()
        palette = [next(cc)['color'] for _ in range(length)]

    else:
        if length <= 28:
            palette = palettes.default_26
        elif length <= len(palettes.default_64):  # 103 colors
            palette = palettes.default_64
        else:
            palette = ['grey' for _ in range(length)]
            logg.info(
                f'the obs value {value_to_plot!r} has more than 103 categories. Uniform '
                "'grey' color will be used for all categories."
            )

    adata.uns[value_to_plot + '_colors'] = palette[:length]


if __name__ == '__main__':
    adata = Adata('./data')
    adata.read_adata('cache/regress.h5ad')
    adata = adata.get_adata()
    gene_expression = adata
    # normalize each cell to have same count number
    # sc.pp.normalize_per_cell(gene_expression)
    # update datastructure to use normalized data
    gene_expression = gene_expression.X

    latent_dim = 50

    # 3. scScope learning
    if gene_expression.shape[0] >= 100000:
        DI_model = DeepImpute.train(
            gene_expression, latent_dim, T=4, batch_size=512, max_epoch=10)
    else:
        DI_model = DeepImpute.train(
            gene_expression, latent_dim, T=4, batch_size=64, max_epoch=300)

    # 4. latent representations and imputed expressions
    latent_code, imputed_val, _ = DeepImpute.predict(
        gene_expression, DI_model)

    # 5. graph clustering
    if latent_code.shape[0] <= 10000:
        label, _, _ = phenograph.cluster(latent_code)
    else:
        label = DeepImpute.scalable_cluster(latent_code)

    # evaluate
    # ARI = adjusted_rand_score(label, label_ground_truth)
    # print(ARI)
    groups = label
    adata.obs['decoder'] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(np.unique(groups).astype('U')),
    )
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    adata.uns['decoder'] = {}
    adata.uns['decoder']['params'] = {'resolution': None, 'random_state': 0}
    _set_default_colors_for_categorical_obs(adata, 'decoder')
    sc.pl.umap(adata, color=['decoder'], size=15)
    sc.tl.draw_graph(adata)
    sc.pl.draw_graph(adata, color=['decoder'], size=15, legend_fontsize=15)

    X_embedded = TSNE(n_components=2).fit_transform(latent_code)

    # visualization of the subpopulation using tSNE
    plt.figure()
    for i in range(10):
        idx = np.nonzero(label == i)[0]
        plt.scatter(X_embedded[idx, 0], X_embedded[idx, 1])
    plt.show()
