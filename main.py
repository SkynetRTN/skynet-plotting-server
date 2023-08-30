import base64
import os
import sys
import traceback
from os import error

from tempfile import mkdtemp
from shutil import rmtree
import numpy as np
from flask import Flask, json, request, jsonify
from flask_cors import CORS
from werkzeug.datastructures import CombinedMultiDict, MultiDict

from cluster_isochrone import get_iSkip, find_data_in_files
from cluster_pro_scraper import scraper_query_object_local, coordinates_to_dist, scraper_query
from gravity_util import find_strain_model_data, find_frequency_model_data, find_bandpass_range
from gaia_util import gaia_match
from plotligo_trimmed import get_data_from_file, bandpassData
from bestFit import fitToData
import uuid
import time

api = Flask(__name__)
CORS(api, origins='http://localhost:3000')
api.debug = True

# Dictionary to store whitened data, time, and last access timestamp
whitened_data = {}

# test
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
    

# Make a new api route that handles strain data updates
####################################################
####################################################
####################################################
## This needs to be fixed to accomodate server side storage

@api.route("/gravitydata", methods=["OPTIONS"])
def handle_options():
    # Set the CORS headers
    headers = {
        "Access-Control-Allow-Origin": "*",  # Replace "*" with the appropriate origin
        "Access-Control-Allow-Methods": "POST, GET, OPTIONS",
        "Access-Control-Allow-Headers": "Content-Type"
    }
    
    return ("", 200, headers)

@api.route("/gravitydata", methods=["POST"])
def get_gravdata():
    try:
        data = request.get_json()
        mass_ratio = float(request.args['ratioMass'])
        total_mass = float(request.args['totalMass'])
        fband = find_bandpass_range(mass_ratio, total_mass)
        session_id = request.args['sessionID']

        if session_id in whitened_data:
            data = whitened_data[session_id]
            strain_whiten = data['whitenedStrain']
            timeData = data['time']
            last_access_time = data['last_access_time']

            # Update the last access timestamp
            whitened_data[session_id]['last_access_time'] = time.time()

            # Process the whitened data and time as needed
            strain_whiten = np.array(strain_whiten)
            timeData = np.array(timeData)
            data = bandpassData(fband[1], fband[0], strain_whiten, timeData)
            midpoint = np.round(data.shape[0] / 2.0)
            buffer = np.ceil(data.shape[0] * 0.25)
            center_of_data = data[int(midpoint - buffer) : int(midpoint + buffer)]
            return jsonify({'data': center_of_data.tolist()})
        else:
            return jsonify({'error': 'Session ID not found'})
    except Exception as e:
        return jsonify({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})
    

@api.route("/gravity", methods=["GET"])
def get_gravity():
    tb = sys.exc_info()[2]
    try:
        mass_ratio = float(request.args['ratioMass'])
        total_mass = float(request.args['totalMass'])
        data = find_strain_model_data(mass_ratio, total_mass)
        for i in range(len(data)):
            data[i][1] = data[i][1] * 10**18
        return json.dumps({'strain_model': data, 'freq_model': find_frequency_model_data(mass_ratio, total_mass)})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


## Changes to gravfile so that whitened data is stored on the server
@api.route("/gravfile", methods=["POST"])
def upload_process_gravdata():
    tempdir = mkdtemp()
    try:
        # Generate a unique identifier for the user's session or data session
        session_id = str(uuid.uuid4())
        # print('ID : ', session_id)
        file = request.files['file']
        file.save(os.path.join(tempdir, "temp-file.hdf5"))
        strain_whiten, timeData = get_data_from_file(os.path.join(tempdir, "temp-file.hdf5"), whiten_data=1)
        timeZero = np.ceil(timeData[0]) + 16
        
        ## The important data seems to always be within the 2 secends preceeding and proceeding the center of the data
        ## rip and tear
        lowInd = 0
        highInd = 0
        for i in range(len(timeData)):
            timeData[i] = timeData[i] - timeZero;
            if -2.1 < timeData[i] < -2.0:
                lowInd = i
            if 1.9 < timeData[i] < 2.0:
                highInd = i   
        print('Before Length: ', len(strain_whiten))
        strain_whiten = strain_whiten[lowInd:highInd]  
        timeData = timeData[lowInd:highInd]  
        print('After Length: ', len(strain_whiten))    

        # Store the whitened data, time, and last access timestamp in the dictionary using the session ID as the key
        whitened_data[session_id] = {
            'whitenedStrain': strain_whiten.tolist(),
            'time': timeData.tolist(),
            'last_access_time': time.time()
        }

        fband = [35, 400]
        if session_id in whitened_data:
            data = whitened_data[session_id]
            strain_whiten = data['whitenedStrain']
            timeData = data['time']
            last_access_time = data['last_access_time']

            # Update the last access timestamp
            whitened_data[session_id]['last_access_time'] = time.time()

            # Process the whitened data and time as needed
            strain_whiten = np.array(strain_whiten)
            timeData = np.array(timeData)
            data = bandpassData(fband[1], fband[0], strain_whiten, timeData)
            midpoint = np.round(data.shape[0] / 2.0)
            buffer = np.ceil(data.shape[0] * 0.25)
            center_of_data = data[int(midpoint - buffer) : int(midpoint + buffer)]
            for i in range(len(center_of_data)):
                center_of_data[i][1] = center_of_data[i][1] * 10**18
        # Set the session ID as a cookie in the response
        response = jsonify({'dataSet': center_of_data.tolist(), 'sessionID': session_id, 'timeZero': timeZero})
        response.set_cookie('session_id', session_id)

        return response
    except Exception as e:
        return jsonify({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})
    finally:
        rmtree(tempdir, ignore_errors=True)



@api.route("/gravprofile", methods=["POST"])
def get_sepctrogram():
    tempdir = mkdtemp()
    try:
        file = request.files['file']
        file.save(os.path.join(tempdir, "temp-file.hdf5"))
        figure, spec_array = get_data_from_file(os.path.join(tempdir, "temp-file.hdf5"), plot_spectrogram=1)
        xbounds = figure.gca().get_xlim()
        ybounds = figure.gca().get_ylim()
        figure.savefig(os.path.join(tempdir, "specplot.png"))

        with open(os.path.join(tempdir, "specplot.png"), "rb") as image2string:
            encoded_image = base64.b64encode(image2string.read())

        return json.dumps({'image': str(encoded_image), 'bounds': str(xbounds)+' '+str(ybounds),
                           'spec_array': np.asarray(spec_array).tolist(), 'x0': str(spec_array.x0),
                           'dx': str(spec_array.dx), 'y0' : str(spec_array.y0), 'dy': str(0.5)})
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
        if not result:
            raise error('No Match in Gaia')
        return json.dumps(result)
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


@api.route("/location-query", methods=["get"])
def get_object_location():
    tb = sys.exc_info()[2]
    try:
        object = request.args['object']
    except:
        raise error({'error': 'Object input invalid type'})
    try:
        result = scraper_query_object_local(object)
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})
    return json.dumps(result)


@api.route("/vizier-query", methods=["post"])
def get_vizier_photometry():
    tb = sys.exc_info()[2]
    try:
        try:
            ra: float = float(request.args['ra'])
            dec: float = float(request.args['dec'])
            r: float = float(request.args['r'])
            coordinates = coordinates_to_dist(ra, dec, r)
            catalog = request.args['catalog']
            file_key = request.args['keys']
            file_data = request.args['data']
            constrain = request.args['constrain']
            if not catalog:
                raise error({'error': 'no catalog!'})
        except Exception as e:
            raise error({'error': 'Object input invalid type'})
        return json.dumps(
            scraper_query(coordinates, constrain, catalog, file_key, file_data)
        ).replace("NaN", "null")
    except Exception as e:
        return json.dumps({'failure': str(e), 'log': traceback.format_tb(e.__traceback__)})


# Periodic cleanup task to remove expired or unused entries
def cleanup_whitened_data():
    expiration_time = 3600 * 3  # Time in seconds after which an entry is considered expired

    current_time = time.time()
    expired_entries = []

    for session_id, data in whitened_data.items():
        last_access_time = data['last_access_time']
        if current_time - last_access_time > expiration_time:
            expired_entries.append(session_id)

    # Remove the expired entries from the whitened_data dictionary
    for session_id in expired_entries:
        del whitened_data[session_id]


def main():
    api.run()

if __name__ == "__main__":
    main()
