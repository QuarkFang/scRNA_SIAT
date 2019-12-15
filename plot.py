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


def show_possibility(list_all: list, names: np.ndarray, scores: np.ndarray, col=4):
    plot_number = len(list_all)
    group_number = np.shape(names)[0]
    plot_row = round(plot_number / col + 0.5)

    list_all = compute_weight(list_all, scores, group_number)

    plt.figure(figsize=(15, 20), dpi=600)

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
        for j in range(group_number):
            plt.bar(x, y[j], bottom=bottom, color=cmap(int(256/(group_number-1))*(j+2)), label=names[i, j])
            bottom = bottom + y[j]

        plt.title('%d vs. rest' % i, size=40)
        plt.legend(loc='upper right', prop={'size': 15})
        plt.xticks(size=15, rotation=90, horizontalalignment='right')
        plt.yticks(size=40)
        plt.grid(axis='y', color='gray', linestyle=':', linewidth=2)
        plt.grid(axis='x', color='gray', linestyle=':', linewidth=2)

        # plt.show()
        plt.tight_layout()
        plt.savefig('groups/group%d.png' % i)
