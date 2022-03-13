from distutils.log import error
import os
from flask import Flask, json, request, render_template
from flask_cors import CORS

# from flask_cors import CORS
from werkzeug.utils import secure_filename
import numpy as np
import ast

from gaia import gaia_args_verify

api = Flask(__name__)
CORS(api)

api.debug = True

cols = [
    "junk",
    "junk",
    "junk",
    "U",
    "B",
    "V",
    "R",
    "I",
    "Jx",
    "Hx",
    "Kx",
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
        return {"error": "Requested data not found"}
    # format data
    try:

        def get_col(number: int):
            return data[:, cols.index(filters[number])]

        r_data = list(
            zip([round(a - b, 4)
                for a, b in zip(get_col(0), get_col(1))], get_col(2))
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


@api.route("/gaia", methods=["POST"])
def get_gaia():
    try:
        data = json.loads(request.get_data())['data']
        result = []
        for star in data:
            result.append(
                {'id': star['id'], 'range': 6.9, 'pm': {'ra': 6.9, 'dec': 6.9}})
    except:
        return json.dumps({"error": "Input invalid type"})

    return json.dumps(result)


def main():
    api.run()


if __name__ == "__main__":
    main()
