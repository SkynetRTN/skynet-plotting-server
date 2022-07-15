import os
import sqlite3
import sys
import traceback
from os import error

from tempfile import mkdtemp
from shutil import rmtree
import numpy as np
from flask import Flask, json, request
from werkzeug.datastructures import CombinedMultiDict, MultiDict
from gravity_util import find_gravity_data
from gaia_util import gaia_match
from plotligo_trimmed import perform_whitening_on_file
from bestFit import fitToData


api = Flask(__name__)

# CORS(api)
# api.debug = True

#test
@api.before_request
def resolve_request_body() -> None:
    ds = [request.args, request.form]
    body = request.get_json()
    if body:
        ds.append(MultiDict(body.items()))

    request.args = CombinedMultiDict(ds)

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


@api.route("/isochrone", methods=["GET"])
def get_data():
    tb = sys.exc_info()[2]
    try:
        age = float(request.args['age'])
        metallicity = float(request.args['metallicity'])
        filters = json.loads(request.args['filters'])
        iSkip = get_iSkip(age, metallicity)
        return json.dumps({'data': find_data_in_files(age, metallicity, filters), 'iSkip': iSkip})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


@api.route("/gravity", methods=["GET"])
def get_gravity():
    tb = sys.exc_info()[2]
    try:
        mass_ratio = float(request.args['ratioMass'])
        total_mass = float(request.args['totalMass'])
        return json.dumps({'data': find_gravity_data(mass_ratio, total_mass)})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


@api.route("/gravfile", methods=["POST"])
def whiten_gravdata():
    # upload_folder = 'temp-grav-data'
    tempdir = mkdtemp()
    try:
        file = request.files['file']
        file.save(os.path.join(tempdir, "temp-file.hdf5"))
        data = perform_whitening_on_file(os.path.join(tempdir, "temp-file.hdf5"))
        midpoint = np.round(data.shape[0]/2.0)
        buffer = np.ceil(data.shape[0] * 0.05)
        center_of_data = data[int(midpoint-buffer): int(midpoint+buffer)]
        return json.dumps({'data': center_of_data.tolist()})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})
    finally:
        rmtree(tempdir, ignore_errors=True)



@api.route("/transient", methods=["POST"])
def get_transient_bestfit():
    tb = sys.exc_info()[2]
    try:
        xdata = request.args['xdata']
        ydata = request.args['ydata']
        filters = request.args['filters']
        guess = request.args['params']
        popt = fitToData(xdata, ydata, filters, guess)
        return json.dumps({'popt': list(popt)})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


@api.route("/gaia", methods=["POST"])
def get_gaia():
    tb = sys.exc_info()[2]
    try:
        try:
            data = request.args['data']
            range = request.args['range']
        except:
            raise error({'error': 'GAIA Input invalid type'})
        result = gaia_match(data, range)
        if result == []:
            raise error('No Match in Gaia')
        return json.dumps(result)
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


def main():
    api.run()


if __name__ == "__main__":
    main()
