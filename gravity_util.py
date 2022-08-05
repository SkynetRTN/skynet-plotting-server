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
