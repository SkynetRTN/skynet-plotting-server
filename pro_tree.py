from sklearn.neighbors import KDTree


def tree_matching(npArray, querys):
    npArray = KDTree(npArray, leaf_size=7000)
    result = []
    for query in querys:
        dist, ind = npArray.query([query], k=1)
        result.append(int(ind[0][0]))
    return result
