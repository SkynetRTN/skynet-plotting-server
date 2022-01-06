from flask import Flask, json, request
import numpy as np
import ast

api = Flask(__name__)


cols = [
    "junk",
    "junk",
    "junk",
    "U",
    "B",
    "V",
    "R",
    "I",
    "J",
    "H",
    "K",
    "uprime",
    "gprime",
    "rprime",
    "iprime",
    "zprime",
    "J",  # TODO: make these distinct
    "H",
    "Ks",
    "junk",
]


def find_data_in_files(
    age: float, metallicity: float, filters: list[str]
) -> list[float]:

    # attempt to retrieve data from files
    try:
        data = np.load(f"./iso-npy-data/Girardi_{age:.2f}_{metallicity:.2f}.npy")
    except FileNotFoundError:
        return {"error": "Requested data not found"}
    # format data
    try:

        def get_col(number: int):
            return data[:, cols.index(filters[number])]

        r_data = list(
            zip([round(a - b, 4) for a, b in zip(get_col(0), get_col(1))], get_col(2))
        )

        # r_data = list(zip(*[data[:, cols.index(i)] for i in filters]))
    except ValueError:
        return {"error": "Requested filter not found"}
    # return
    return r_data


@api.route("/isochrone", methods=["GET"])
def get_data():
    try:
        age = float(request.args.get("age"))
        metallicity = float(request.args.get("metallicity"))
        filters = ast.literal_eval(request.args.get("filters"))
        print(filters)
    except ValueError:
        return json.dumps({"error": "Input invalid type"})
    return json.dumps(find_data_in_files(age, metallicity, filters))


def main():
    api.run()


if __name__ == "__main__":
    main()
