from sklearn.neighbors import KDTree
from scipy.spatial import KDTree as KDTree_sci


def tree_matching(npArray, querys):
    tree = KDTree(npArray, leaf_size=7000)
    result = []
    for query in querys:
        dist, ind = tree.query([query], k=1)
        result.append(int(ind[0][0]))
    return result


def tree_matching_sci(npArray, querys):
    tree = KDTree_sci(npArray, 2)
    return tree.query(querys)[1]
