import matplotlib.pyplot as plt
import numpy as np


def compute_weight(list_all: list, scores: np.ndarray, group_number: int):
    for i, list_row in enumerate(list_all):
        for j, list_col in enumerate(list_row):
            plot_weight = []
            for k in range(group_number):
                if k in list_col[1]:
                    plot_weight.append(scores[i, k])
                else:
                    plot_weight.append(0.0)
            list_col.append(plot_weight)
    return list_all


def show_possibility(list_all: list, names: np.ndarray, scores: np.ndarray, gene_number=10, col=4):
    plot_number = len(list_all)
    group_number = np.shape(names)[0]
    plot_row = round(plot_number / col + 0.5)

    list_all = compute_weight(list_all, scores, gene_number)

    plt.figure(figsize=(15, 20), dpi=200)

    cmap = plt.cm.get_cmap("Spectral")

    for i, list_row in enumerate(list_all):
        plt.cla()
        x = []
        y = []
        for j, list_col in enumerate(list_row):
            x.append(list_col[0])
            y.append(list_col[3])
        y = np.array(y).T
        bottom = np.zeros(len(list_row))
        for j in range(gene_number):
            plt.bar(x, y[j], width=0.6, bottom=bottom, color=cmap(int(256/(gene_number-1))*j), label=names[i, j])
            bottom = bottom + y[j]

        plt.title('Group %d' % i, size=40)
        plt.legend(loc='upper right', prop={'size': 15})
        plt.xticks(size=15, rotation=90, horizontalalignment='right')
        plt.yticks(size=20)
        plt.grid(axis='y', color='gray', linestyle=':', linewidth=2)
        plt.grid(axis='x', color='gray', linestyle=':', linewidth=2)

        # plt.show()
        plt.tight_layout()
        plt.savefig('groups/group%d.png' % i)
