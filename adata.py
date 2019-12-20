import numpy as np
import scanpy as sc
import os
import anndata
from natsort import natsorted
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import phenograph
import scscope as DeepImpute
import pandas as pd


class Adata:
    dpi = 200
    save_dpi = 200
    adata = None
    result_file = './write/pbmc3k.h5ad'

    def __init__(self, filepath: str):
        sc.settings.verbosity = 3
        sc.logging.print_versions()
        sc.settings.set_figure_params(dpi=self.dpi, dpi_save=self.save_dpi)

        self.adata = sc.read_10x_mtx(filepath, var_names='gene_symbols', cache=True)
        self.adata.var_names_make_unique()

    def get_adata(self):
        return self.adata

    def set_adata(self, adata: anndata.AnnData):
        self.adata = adata

    def set_dpi(self, dpi=200, save_dpi=200):
        self.dpi = dpi
        self.save_dpi = save_dpi
        sc.settings.set_figure_params(dpi=dpi, dpi_save=save_dpi)

    def set_categories(self, cat: list):
        self.adata.obs['louvain'].cat.categories = cat

    def save_adata(self, filename):
        path = os.fspath(filename)
        self.adata.write(path)

    def read_adata(self, filename):
        self.adata = sc.read_h5ad(filename)

    def filter(self, min_genes=200, min_cells=3):
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)

    def mit(self, max_genes=4000, max_percent=0.05):
        mito_genes = self.adata.var_names.str.startswith('MT-')

        self.adata.obs['percent_mito'] = np.sum(
            self.adata[:, mito_genes].X, axis=1).A1 / np.sum(self.adata.X, axis=1).A1

        self.adata.obs['n_counts'] = self.adata.X.sum(axis=1).A1

        self.adata = self.adata[self.adata.obs['n_genes'] < max_genes, :]
        self.adata = self.adata[self.adata.obs['percent_mito'] < max_percent, :]

    def normalize(self):
        sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=1e4)

        sc.pp.log1p(self.adata)

        self.adata.raw = self.adata

    def regress(self):
        sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
        # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        # sc.pl.highly_variable_genes(adata)

        self.adata = self.adata[:, self.adata.var['highly_variable']]

        sc.pp.regress_out(self.adata, ['n_counts', 'percent_mito'], n_jobs=None)
        sc.pp.scale(self.adata, max_value=10)

    def pca(self, show=False):
        sc.tl.pca(self.adata, svd_solver='arpack')
        if show:
            sc.pl.pca_variance_ratio(self.adata, log=True)

    def neighbors(self, n_neighbors=5, n_pcs=20):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    def umap(self, show=False, save=False, color=None):
        if color is None:
            color = ['louvain']
        sc.tl.umap(self.adata)
        sc.tl.louvain(self.adata)
        if show:
            sc.pl.umap(self.adata, color=color, size=15)
        if save and color == ['louvain']:
            sc.pl.umap(self.adata, color=color, size=15, save='_groups.png')
        elif save:
            sc.pl.umap(self.adata, color=color, size=15, save='_groups_genes.png')

    def rank_genes_groups(self, n_genes=25, show=False, save=False):
        sc.tl.rank_genes_groups(self.adata, 'louvain', method='wilcoxon', log_transformed=True)
        np.save('adata_scores', self.adata.uns['rank_genes_groups']['scores'])
        np.save('adata_names', self.adata.uns['rank_genes_groups']['names'])
        if show:
            sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False)
        if save:
            sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False, save='_top25.png')

    def regroup(self, index: list, show=False, save=False):
        """
        Parameters
        ----------

        index: 'list'
            like [(1, 2), (3, 6, 9)]
        show: 'bool'
        save: 'bool'
        """

        group_number = np.shape(self.adata.uns['louvain_colors'])[0]
        index_dict = {}
        pointer = 0
        for t in index:
            for group in t:
                index_dict[group] = pointer
            pointer = pointer + 1
        for i in range(group_number):
            if i not in index_dict:
                index_dict[i] = pointer
                pointer = pointer + 1

        louvain = self.adata.obs['louvain']
        for i, _ in enumerate(louvain):
            louvain[i] = str(index_dict[int(louvain[i])])

        self.adata.uns['louvain_colors'] = self.adata.uns['louvain_colors'][0:pointer]
        self.adata.obs['louvain'].cat.categories = self.adata.obs['louvain'].cat.categories[0:pointer]

        if show:
            sc.pl.umap(self.adata, color=['louvain'], size=15)
        if save:
            sc.pl.umap(self.adata, color=['louvain'], size=15, save='_regroups.png')

    def paga(self, show=False, save=False, color=None):
        if color is None:
            color = ['louvain']
        sc.tl.paga(self.adata, groups='louvain')
        if show:
            sc.pl.paga(self.adata, color=color, title="", fontsize=10)
        if save and color == ['louvain']:
            sc.pl.paga(self.adata, color=color, title="", fontsize=10, save='_groups.png')
        elif save:
            sc.pl.paga(self.adata, color=color, title="", fontsize=10, save='_groups_genes.png')

    def draw_graph(self, show=False, save=False, color=None):
        if color is None:
            color = ['louvain']
        sc.tl.draw_graph(self.adata)
        if show:
            sc.pl.draw_graph(self.adata, color=color, size=15, legend_fontsize=15)
        if save and color == ['louvain']:
            sc.pl.draw_graph(self.adata, color=color, size=15, save='_groups.png', legend_fontsize=15)
        elif save:
            sc.pl.draw_graph(self.adata, color=color, size=15, save='_groups_genes.png', legend_fontsize=15)

    def _set_default_colors_for_categorical_obs(self, value_to_plot):
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

        categories = self.adata.obs[value_to_plot].cat.categories
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

        self.adata.uns[value_to_plot + '_colors'] = palette[:length]

    def dp_group(self, show=False, show_tsne=False, show_fr=False):
        """
        using deep learning method to classify

        """
        gene_expression = self.adata.X
        # normalize each cell to have same count number
        # sc.pp.normalize_per_cell(gene_expression)
        # update data structure to use normalized data
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
        self.adata.obs['decoder'] = pd.Categorical(
            values=groups.astype('U'),
            categories=natsorted(np.unique(groups).astype('U')),
        )
        sc.tl.pca(self.adata, svd_solver='arpack')
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)
        self.adata.uns['decoder'] = {}
        self.adata.uns['decoder']['params'] = {'resolution': None, 'random_state': 0}
        self._set_default_colors_for_categorical_obs('decoder')

        if show:
            sc.pl.umap(self.adata, color=['decoder'], size=15)
        if show_tsne:
            X_embedded = TSNE(n_components=2).fit_transform(latent_code)

            # visualization of the subpopulation using tSNE
            plt.figure()
            for i in range(10):
                idx = np.nonzero(label == i)[0]
                plt.scatter(X_embedded[idx, 0], X_embedded[idx, 1])
            plt.show()
        if show_fr:
            sc.tl.draw_graph(self.adata)
            sc.pl.draw_graph(self.adata, color=['decoder'], size=15, legend_fontsize=15)