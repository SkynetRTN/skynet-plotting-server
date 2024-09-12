import base64
import os
import sys
import redis
import pickle
import traceback
from os import error
import pandas as pd
import gwpy.timeseries
import gwpy.frequencyseries
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.filter import resample_to_delta_t, highpass
from pycbc.types.frequencyseries import FrequencySeries
from pycbc.types.timeseries import TimeSeries
from tempfile import mkdtemp
from shutil import rmtree
from pycbc.filter import matched_filter
import numpy as np
from flask import Flask, json, request, jsonify
from werkzeug.datastructures import CombinedMultiDict, MultiDict

from cluster_isochrone import get_iSkip, find_data_in_files
from cluster_pro_scraper import scraper_query_object_local, coordinates_to_dist, scraper_query
from gravity_util import find_strain_model_data, find_frequency_model_data, find_bandpass_range, find_normalization, \
   find_raw_fmodel
from gaia_util import gaia_match
from plotligo_trimmed import get_data_from_file, bandpassData
from transient_model_fit import fit
import uuid
import time
from flask_cors import CORS

api = Flask(__name__)
api.debug = True
CORS(api)

r = redis.StrictRedis(host='localhost', port=6379, db=0)
DATA_EXPIRATION = 86400

# Dictionary to store whitened data, time, and last access timestamp, and raw data (maybe)
#whitened_data = {}


# test
@api.before_request
def resolve_request_body() -> None:
    ds = [request.args, request.form]
    body = request.get_json(force=True, silent=True)
    if body:
        ds.append(MultiDict(body.items()))

    request.args = CombinedMultiDict(ds)

    # ds = {**request.args, **request.form}
    # body = request.get_json(force=True, silent=True)
    # if body:
    #     ds.append(**MultiDict(body.items()))

    # request.args = CombinedMultiDict(ds)

#CLUSTER
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
        "Access-Control-Allow-Headers": "Content-Type",
    }

    return ("", 200, headers)

#Pulls model based on user selected params, processes data and sends it back. (purple wave)
@api.route("/gravitydata", methods=["POST"])
def get_gravdata():
    try:
        data = request.get_json()
        mass_ratio = float(request.args['ratioMass'])
        total_mass = float(request.args['totalMass'])
        phase = float(request.args['phaseMass'])
        fband = find_bandpass_range(mass_ratio, total_mass, phase)
        hp = find_raw_fmodel(mass_ratio, total_mass, phase)
        session_id = request.args['sessionID']
        
        if r.exists(session_id):
            data = pickle.loads(r.get(session_id))
            strain_whiten = data['whitenedStrain']
            timeData = data['time']
            last_access_time = data['last_access_time']
            rawTimeseries = data['rawTimeseries']
            rawTimeseries = resample_to_delta_t(highpass(rawTimeseries, fband[0]), 1.0 / 2048)
            conditioned = rawTimeseries.crop(2, 2)
            psd = conditioned.psd(4)
            psd = interpolate(psd, conditioned.delta_f)
            psd = inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),
                                              low_frequency_cutoff=fband[0])
            ## Now that we have the psd, rawtimeseries, and a model to compare them to
            ## we can use the pycbc filter matching and find an approrpiate SNR for all of these!!
            hp = TimeSeries(hp.real, delta_t=1 / 16384)
            hp = resample_to_delta_t(hp, 1.0 / 2048)
            hp.resize(len(conditioned))
            template = hp.cyclic_time_shift(hp.start_time)
            snr = matched_filter(template,
                                 conditioned,
                                 psd=psd,
                                 low_frequency_cutoff=20)
            ## may not need as a timeseries -- just find max SNR I guess
            snr = snr.crop(4 + 4, 4)
            peak = abs(snr).numpy().argmax()
            snrMax = abs(snr[peak])
            print('SNR Max is: ', snrMax)
            # Update the last access timestamp
            # Process the whitened data and time as needed
            ## there is a VERY strong chance that we should just be storing the whitened strain data in fequency space to save one transform
            ## pronorm = find_normalization(mass_ratio, total_mass)
            strain_whiten = np.array(strain_whiten)

            # Now I think we can maybe just apply it here

            for i in range(len(strain_whiten)):
                strain_whiten[i] = strain_whiten[i] * snrMax



            ## Nt = len(strain_whiten)
            ## strain_whiten = np.fft.rfft(strain_whiten)
            ## strain_whiten = strain_whiten * pronorm
            ## strain_whiten = np.fft.irfft(strain_whiten, n=Nt)
            timeData = np.array(timeData)
            data = bandpassData(fband[1], fband[0], strain_whiten, timeData)
            midpoint = np.round(data.shape[0] / 2.0)
            buffer = np.ceil(data.shape[0] * 0.25)
            center_of_data = data[int(midpoint - buffer): int(midpoint + buffer)]

            if np.isnan(center_of_data[0][1]) == True:
                center_of_data = np.nan_to_num(center_of_data, nan=0.0)

            r.expire(session_id, DATA_EXPIRATION)
            print('Data trouble: ', center_of_data[0][1])

            # for i in range(len(center_of_data)):
            #     center_of_data[i][1] = center_of_data[i][1]
            return jsonify({'data': center_of_data.tolist(), 'snrMax': snrMax})
        else:
            return jsonify({'error': 'Session ID not found'})
    except Exception as e:
        return jsonify({'err': str(e), 'log': traceback.format_tb(e.__traceback__)})


# takes user defined params and returns a model (orange wave and yellow bounds)
@api.route("/gravity", methods=["GET"])
def get_gravity():
    tb = sys.exc_info()[2]
    try:
        mass_ratio = float(request.args['ratioMass'])
        total_mass = float(request.args['totalMass'])
        mass_ratioStrain = float(request.args['ratioMassStrain'])
        total_massStrain = float(request.args['totalMassStrain'])
        phase = float(request.args['phaseStrain'])
        data = find_strain_model_data(mass_ratioStrain, total_massStrain, phase)
        for i in range(len(data)):
            data[i][1] = data[i][1] * 10 ** 22.805
            if i > 0.65 * len(data):
                data[i][1] = 0

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
        print('ID : ', session_id)
        if request.args['default_set'] == 'true':
            strain_whiten, timeData, PSD, rawTimeseries, timeOfRecord, timeZero = get_data_from_file(
            './gravity-model-data/default-set/GW150914.hdf5'
                , whiten_data=1)    
        else:
            file = request.files['file']
            file.save(os.path.join(tempdir, "temp-file.hdf5"))

            strain_whiten, timeData, PSD, rawTimeseries, timeOfRecord, timeZero = get_data_from_file(
                os.path.join(tempdir, "temp-file.hdf5"), whiten_data=1)
        ## We need to export all the string info from the file upon loading for naming purposes

        print('Time Data Raw: ', timeData)
        print('Raw TimeSeries Data is: ', rawTimeseries)
        rawTimeseries = TimeSeries(rawTimeseries.to_pycbc(), delta_t=1 / 16384)
        ## The important data seems to always be within the 2 secends preceeding and proceeding the center of the data
        ## rip and tear
        lowInd = 0
        highInd = 0
        for i in range(len(timeData)):
            timeData[i] = timeData[i] - timeZero;
            if 13.9 < timeData[i] < 14.0:
                lowInd = i
            if 17.9 < timeData[i] < 18.0:
                highInd = i
        print('Before Length: ', len(strain_whiten))
        strain_whiten = strain_whiten[lowInd:highInd]
        timeData = timeData[lowInd:highInd]
        print('TimeData cropped: ', timeData)
        print('After Length: ', len(strain_whiten))

        # Store the whitened data, time, and last access timestamp in the dictionary using the session ID as the key
        # for SNR reasons we may also store the raw data timeseries and its PSD here too
        data  = {
            'whitenedStrain': strain_whiten.tolist(),
            'time': timeData.tolist(),
            'last_access_time': time.time(),
            'rawTimeseries': rawTimeseries,
            'PSD': PSD
        }
        r.set(session_id, pickle.dumps(data))
        r.expire(session_id, DATA_EXPIRATION)

        fband = [35, 400]
        if True:
            strain_whiten = data['whitenedStrain']
            timeData = data['time']
            last_access_time = data['last_access_time']


            # Process the whitened data and time as needed
            strain_whiten = np.array(strain_whiten)
            timeData = np.array(timeData)
            data = bandpassData(fband[1], fband[0], strain_whiten, timeData)
            midpoint = np.round(data.shape[0] / 2.0)
            buffer = np.ceil(data.shape[0] * 0.25)
            center_of_data = data[int(midpoint - buffer): int(midpoint + buffer)]
            # for i in range(len(center_of_data)):
            #     center_of_data[i][1] = center_of_data[i][1]
        # Set the session ID as a cookie in the response
        response = jsonify({'dataSet': center_of_data.tolist(), 'sessionID': session_id, 'timeZero': timeZero.tolist(),
                            'timeOfRecord': timeOfRecord})
        response.set_cookie('session_id', session_id)

        return response
    except Exception as e:
        print(e)
        return jsonify({ 'err': str(e), 'log': traceback.format_tb(e.__traceback__)}), 400
    finally:
        rmtree(tempdir, ignore_errors=True)


@api.route("/gravprofile", methods=["POST"])
def get_sepctrogram():
    tempdir = mkdtemp()
    try:
        if request.args['default_set'] == 'true':
            figure, spec_array = get_data_from_file(
            './gravity-model-data/default-set/GW150914.hdf5'
                , plot_spectrogram=1)
        else:
            file = request.files['file']
            file.save(os.path.join(tempdir, "temp-file.hdf5"))
            figure, spec_array = get_data_from_file(os.path.join(tempdir, "temp-file.hdf5"), plot_spectrogram=1)
        # NetworkSNR = NetworkSNR.max()
        # print('The Network SNR for this one is: ', NetworkSNR)
        xbounds = figure.gca().get_xlim()
        ybounds = figure.gca().get_ylim()
        figure.savefig(os.path.join(tempdir, "specplot.png"))

        with open(os.path.join(tempdir, "specplot.png"), "rb") as image2string:
            encoded_image = base64.b64encode(image2string.read())

        print(sys.getdefaultencoding())

        with open("zServerImage.png", 'wb') as f:
            img = base64.b64decode(bytes(str(encoded_image, encoding=sys.getdefaultencoding()), encoding=sys.getdefaultencoding()))
            f.write(img)

        return json.dumps({'image': str(encoded_image), 'bounds': str(xbounds) + ' ' + str(ybounds),
                           'spec_array': np.asarray(spec_array).tolist(), 'x0': str(spec_array.x0),
                           'dx': str(spec_array.dx), 'y0': str(spec_array.y0), 'dy': str(0.5)})
    except Exception as e:
        return json.dumps({'err': str(e), 'log': traceback.format_tb(e.__traceback__)}), 400
    finally:
        rmtree(tempdir, ignore_errors=True)

# NON GRAVITY ENDPOINTS

@api.route("/transient", methods=["POST"])
def get_transient_bestfit():
    tb = sys.exc_info()[2]
    try:
        popt = fit(request.args['xdata'], request.args['ydata'], request.args['filters'], request.args['params'])
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


# # Periodic cleanup task to remove expired or unused entries
# def cleanup_whitened_data():
#     expiration_time = 3600 * 3  # Time in seconds after which an entry is considered expired

#     current_time = time.time()
#     expired_entries = []

#     for session_id, data in whitened_data.items():
#         last_access_time = data['last_access_time']
#         if current_time - last_access_time > expiration_time:
#             expired_entries.append(session_id)

#     # Remove the expired entries from the whitened_data dictionary
#     for session_id in expired_entries:
#         del whitened_data[session_id]


def main():
    # Run the cleanup task every hour (3600 seconds)
    # api.before_first_request(cleanup_whitened_data)
    api.run(port=5001)


if __name__ == "__main__":
    main()
