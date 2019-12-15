import dataset as db
import match
import numpy as np
import plot as plt


if __name__ == '__main__':
    # dataset.add(6, 'Human', 'Kidney', 'UBERON_0002113', 'Normal', 'Normal cell',
    #             'Proximal tubular cell', 'NA', 'AA1',
    #             '9263997')

    # a = db.query('%NFKBIA%')
    # print(a)

    b = np.load('adata_scores.npy')
    a = np.load('adata_names.npy')
    c, names, scores = match.match(a, b, top_gene=25)
    plt.show_possibility(c, names, scores, gene_number=25)  # gene_number < top_gene
