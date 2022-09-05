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
                f"freqmodel-{total_mass:.3f}-{mass_ratio:.3f}.dat"))
    model[:,0] = merger_time + model[:,0]
    # for cord in model:
    #     print(extract_value_in_spec(cord, spectrogram))
    return model

def extract_value_in_spec(x, spect):
    return spect[x[0], x[1]]

# spec = np.asarray([[1,2,3], [4, 5, 6], [7, 8, 9]])
#
# print(extract_model_from_spectrogram(1.000, 9.953, 5 ,spec)[1:])