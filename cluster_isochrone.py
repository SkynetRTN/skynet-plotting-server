import os
import sqlite3
from os import error
import numpy as np

cols = [
    "junk",
    "junk",
    "U",
    "B",
    "V",
    "R",
    "I",
    "uprime",
    "gprime",
    "rprime",
    "iprime",
    "zprime",
    "J",
    "H",
    "K",
    "W1",
    "W2",
    "W3",
    "W4",
    "G",
    "BP",
    "RP",
    "junk",
]


def find_data_in_files(age: float, metallicity: float, filters: list) -> list:

    # attempt to retrieve data from files
    try:
        data = np.load(
            os.path.join(
                os.path.dirname(__file__),
                "iso-npy-data",
                f"Girardi_{age:.2f}_{metallicity:.2f}.npy",
            )
        )

    except FileNotFoundError:
        raise ValueError({"error": "Requested data not found"})
    # format data
    try:

        def get_col(number: int):
            return data[:, cols.index(filters[number])]

        r_data = list(
            zip([round(a - b, 4)
                 for a, b in zip(get_col(0), get_col(1))], get_col(2))
        )

        # r_data = list(zip(*[data[:, cols.index(i)] for i in filters]))
    except:
        raise error({"error": "Requested filter not found"})
    # return
    return r_data


def get_iSkip(age, metallicity):
    iSkip_file = os.path.join(os.path.dirname(__file__), 'iSkip.sqlite')
    conn = sqlite3.connect(iSkip_file)
    result = -1
    try:
        cur = conn.cursor()
        sources = cur.execute(
            'SELECT * FROM iSkip_forsql WHERE age = ? AND metallicity = ?', [age, metallicity]).fetchall()
        if sources != []:
            result = sources[0][2]
    except:
        raise error({"error": "Cannot Pull from iSkip Database"})
    finally:
        conn.close()
    return result
