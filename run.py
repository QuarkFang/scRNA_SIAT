import dataset as db
import match
import numpy as np
import plot as plt
from adata import Adata


if __name__ == '__main__':

    adata = Adata('./data')
    # adata.set_dpi(save_dpi=400)
    # adata.filter()
    # adata.mit()
    # adata.normalize()
    # adata.regress()
    # adata.pca()
    # adata.neighbors()
    # adata.umap(save=True)
    # adata.umap(color=['KRT8', 'KRT18'], save=True)
    # adata.save_adata('cache/umap.h5ad')
    adata.read_adata('cache/umap.h5ad')
    adata.rank_genes_groups()
    b = np.load('adata_scores.npy')
    a = np.load('adata_names.npy')
    c, names, scores = match.match(a, b, top_gene=25)
    plt.show_possibility(c, names, scores, gene_number=25)  # gene_number < top_gene
    adata.regroup([(0, 1, 3, 4), (2, 6, 7, 9, 10, 13, 14)], save=True)
    adata.set_categories(['0/Monocyte', '1/NKT/B', '2/CD1C-CD141',
                          '3/Unknown', '4/Unknown', '5/Unknown',
                          '6/Plasma'])
    adata.draw_graph(save=True)
    adata.draw_graph(color=['CD79A', 'CD3D'], save=True)
    adata.paga()
