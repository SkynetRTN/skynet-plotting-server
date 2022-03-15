import time
import numpy as np
from sklearn.neighbors import KDTree

duration = 0
build = 0
for _ in range(900):
    rng = np.random.RandomState(0)
    X = rng.random_sample((30000, 2))  # 10 points in 3 dimensions
    x = rng.random_sample((1, 2))
    start = time.time()
    tree = KDTree(X, leaf_size=2)
    end = time.time()
    build += (end - start)
    start = time.time()
    dist, ind = tree.query(x, k=1)
    end = time.time()
    duration += (end - start)
# print(ind)  # indices of 3 closest neighbors
print(duration)
print(build)
