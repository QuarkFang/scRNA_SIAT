import dataset as db
import numpy as np
from operator import itemgetter


def recarray2ndarray(rec: np.recarray):
    res = []
    for row in rec:
        res.append(list(row))
    return np.array(res).T


def marker_cell(names: np.ndarray, scores: np.ndarray, top_gene=10):
    names = names[:, 0:top_gene]
    scores = scores[:, 0:top_gene]
    group_num = np.shape(names)[0]

    all_names = []
    for row in names:
        cell_row = []
        for marker in row:
            query = db.query(marker)
            cell_name = []
            for cell in query:
                cell_name.append(cell['species_type'] + ',\n' + cell['tissue_type'] + ',\n' + cell['cell_name'])
            cell_row.append(cell_name)
        all_names.append(cell_row)

    all_union_names = []
    for row in all_names:
        union_names = []
        for cell_list in row:
            union_names = list(set(union_names).union(set(cell_list)))
        all_union_names.append(union_names)

    list_a = []  # Name
    list_b = []  # Name->Marker index
    list_c = []  # Name->Score
    for i, row in enumerate(all_union_names):
        list_a_row = []
        list_b_row = []
        list_c_row = []
        for cell_name in row:
            list_a_row.append(cell_name)

            score = 0
            list_b_col = []
            for j, col in enumerate(all_names[i]):
                if cell_name in col:
                    list_b_col.append(j)
                    score = score + scores[i, j]
            list_b_row.append(list_b_col)
            list_c_row.append(score)

        list_a.append(np.array(list_a_row))
        list_b.append(np.array(list_b_row))
        list_c.append(np.array(list_c_row))
    return list_a, list_b, list_c


def cmp(x, y):
    return x > y


def list_sort(list_a: list, list_b: list, list_c: list):
    list_all = []
    for i, row in enumerate(list_a):
        list_all_row = []
        for j, col in enumerate(row):
            list_all_col = [list_a[i][j], list_b[i][j], list_c[i][j]]
            list_all_row.append(list_all_col)
        list_all_row.sort(key=itemgetter(2), reverse=True)

        list_len = len(list_all_row)
        if list_len >= 10:
            list_all_row = list_all_row[0:10]
        else:
            for k in range(10-list_len):
                list_all_row.append(['NA%d,\nto supplement the form' % (k+1), [], 0.0])
        list_all.append(list_all_row)

    return list_all


def normalize(list_all: list):
    for list_row in list_all:
        sum_weight = 0
        for list_col in list_row:
            sum_weight = sum_weight + list_col[2]

        for list_col in list_row:
            list_col[2] = list_col[2] / sum_weight

    return list_all


def match(names: np.recarray, scores: np.recarray, top_gene=25):
    names = recarray2ndarray(names)
    scores = recarray2ndarray(scores)
    list_a, list_b, list_c = marker_cell(names, scores, top_gene=top_gene)
    list_all = list_sort(list_a, list_b, list_c)
    return list_all, names, scores
