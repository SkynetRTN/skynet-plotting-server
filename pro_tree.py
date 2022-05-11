from sklearn.neighbors import KDTree
from scipy.spatial import KDTree as KDTree_sci
import numpy as np
import grispy as gsp


def tree_matching(npArray, querys):
    tree = KDTree(npArray, leaf_size=7000)
    result = []
    for query in querys:
        dist, ind = tree.query([query], k=1)
        result.append(int(ind[0][0]))
    return result


def tree_matching_sci(npArray, querys):
    newArray = []
    for data in npArray:
        newArray.append([data[1], data[0]])
    newQuerys = []
    for query in querys:
        newQuerys.append([query[1], query[0]])
    tree = KDTree_sci(newArray, 2)
    return tree.query(newQuerys)[1]
    tree = KDTree_sci(npArray, 2)
    return tree.query(querys)[1]


def tree_matching_grispy(array, querys):

    grid = gsp.GriSPy(np.array(array), N_cells=128)
    near_dist, near_ind = grid.nearest_neighbors(np.array(querys), n=1)
    return [near_dist, near_ind]
