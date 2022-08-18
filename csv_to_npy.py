import os

import numpy as np
# import pandas as pd

dir = 'iso-npy-data-beta'


for file in sorted(os.listdir(os.path.join(os.path.dirname(__file__), dir))):
    if 'csv' in str(file):
        # print(file)
        arr = np.genfromtxt(os.path.join(os.path.dirname(__file__), dir, file), delimiter=",")
        # print(arr)
        np.save(os.path.join(os.path.dirname(__file__), dir, file.split(".csv")[0]), arr)