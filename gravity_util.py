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
        step = data[1, 0].tolist() - data[0,0].tolist()
        end = data[-1, 0].tolist()
        start = data[0, 0].tolist()
        data = data[:, :2].tolist()
        # we want to add some elements to extend the dataset and make it more stable, I'm not going to start 
        # very accurate, but lets see what we can do
        # for i in range(100000):
        #     data.insert(-1, [end + (step * i), 0])
        #     data. insert(0, [start - (step * i), 0])
        return data

    except FileNotFoundError:
        raise ValueError({"error": "Requested strain model not found"})

# pull the bandpass frequencies from the files

def find_bandpass_range(mass_ratio, total_mass):
    try:
        data = np.genfromtxt(
            os.path.join(
                os.path.dirname(__file__),
                "gravity-model-data", "strain-model", f"bf_{total_mass:0.3f}",
                f"bandpassData-{total_mass:.3f}-{mass_ratio:.3f}.dat"
            ),
            skip_header=1)
        data = data.tolist()
        return data

    except FileNotFoundError:
        raise ValueError({"error": "Requested bandpass data not found"})


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


