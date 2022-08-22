import os
import sqlite3
import sys
import traceback
from os import error

from tempfile import mkdtemp
from shutil import rmtree
import numpy as np
from flask import Flask, json, request, send_file
from werkzeug.datastructures import CombinedMultiDict, MultiDict
from gravity_util import find_strain_model_data, find_frequency_model_data
import ast
from cluster_isochrone import get_iSkip, find_data_in_files, find_data_in_files_beta
from gravity_util import find_gravity_data
from gaia import gaia_args_verify
from gaia_util import gaia_match
from plotligo_trimmed import get_data_from_file
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


@api.route("/isochrone-beta", methods=["GET"])
def get_data_beta():
    tb = sys.exc_info()[2]
    try:
        age = float(request.args['age'])
        metallicity = float(request.args['metallicity'])
        filters = json.loads(request.args['filters'])
        iSkip = get_iSkip(age, metallicity)
        return json.dumps({'data': find_data_in_files_beta(age, metallicity, filters), 'iSkip': -1})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


@api.route("/gravity", methods=["GET"])
def get_gravity():
    tb = sys.exc_info()[2]
    try:
        mass_ratio = float(request.args['ratioMass'])
        total_mass = float(request.args['totalMass'])
        return json.dumps({'strain_model': find_strain_model_data(mass_ratio, total_mass), 'freq_model': find_frequency_model_data(mass_ratio, total_mass)})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


@api.route("/gravfile", methods=["POST"])
def upload_process_gravdata():
    # upload_folder = 'temp-grav-data'
    tempdir = mkdtemp()
    try:
        file = request.files['file']
        file.save(os.path.join(tempdir, "temp-file.hdf5"))
        data = get_data_from_file(os.path.join(tempdir, "temp-file.hdf5"), whiten_data=1)
        midpoint = np.round(data.shape[0]/2.0)
        buffer = np.ceil(data.shape[0] * 0.05)
        center_of_data = data[int(midpoint-buffer): int(midpoint+buffer)]
        return json.dumps({'data': center_of_data.tolist()})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})
    finally:
        rmtree(tempdir, ignore_errors=True)


@api.route("/gravprofile", methods=["POST"])
def get_sepctrogram():
    tempdir = mkdtemp()
    try:
        file = request.files['file']
        file.save(os.path.join(tempdir, "temp-file.hdf5"))
        figure = get_data_from_file(os.path.join(tempdir, "temp-file.hdf5"), plot_spectrogram=1)
        xbounds = figure.gca().get_xlim()
        ybounds = figure.gca().get_ylim()
        figure.savefig(os.path.join(tempdir, "specplot.png"))
        ret = send_file(os.path.join(tempdir, "specplot.png"), mimetype='image/png')
        ret.headers['bounds'] = str(xbounds)+' '+str(ybounds)
        ret.headers['Access-Control-Expose-Headers'] = 'bounds'
        return ret
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
