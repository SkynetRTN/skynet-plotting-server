import numpy as np
import os

def find_strain_model_data(mass_ratio, total_mass):
    try:
        data = np.genfromtxt(
            os.path.join(
                os.path.dirname(__file__),
                "gravity-model-data", "strain-model", f"mt_{total_mass:0.3f}",
                f"gravdata-{total_mass:.3f}-{mass_ratio:.3f}.dat"
            ),
            skip_header=1)
        return data[:,:2].tolist()

    except FileNotFoundError:
        raise ValueError({"error": "Requested strain model not found"})


def find_frequency_model_data(mass_ratio, total_mass):
    try:
        data = np.genfromtxt(
            os.path.join(
                os.path.dirname(__file__),
                "gravity-model-data", "frequency-model", f"mt_{total_mass:0.3f}",
                f"freqmodel-{total_mass:.3f}-{mass_ratio:.3f}.dat"
            ),
            skip_header=1)
        return data[:,:2].tolist()

    except FileNotFoundError:
        raise ValueError({"error": "Requested frequency model not found"})

def extract_model_from_spectrogram(mass_ratio, total_mass, merger_time, spectrogram):
    model = np.genfromtxt(
            os.path.join(
                os.path.dirname(__file__),
                "gravity-model-data", "frequency-model", f"mt_{total_mass:0.3f}",
                f"freqmodel-{total_mass:.3f}-{mass_ratio:.3f}.dat"), skip_header=1)
    model[:,0] = merger_time + model[:,0]
    extracted_output = np.zeros(np.shape(model))
    i = 0
    for vals in model:
        extracted_output[i][0] = vals[0]
        extracted_output[i][1] = extract_value_in_spec(vals,spectrogram)
        i = i+1
    return extracted_output

def extract_value_in_spec(x, spect):
    x0 = 1126259461.4
    # x0 = 0.14
    dx = 0.003199999809265137
    xcord = int((x[0] - x0) / dx)

    y0 = 11.25395
    dy = 0.5
    ycord = int((x[1] - y0) / dy)
    return spect[xcord, ycord]

# spec = np.asarray([[1,2,3], [4, 5, 6], [7, 8, 9]])
#
# print(extract_model_from_spectrogram(1.000, 9.953, 5 ,spec)[1:])