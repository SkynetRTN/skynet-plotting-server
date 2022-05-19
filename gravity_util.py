import numpy as np
import os

def find_gravity_data(mass_ratio, total_mass):
    try:
        data = np.genfromtxt(
            os.path.join(
                os.path.dirname(__file__),
                "gravity-data", f"mt_{total_mass:0.3f}",
                f"gravdata-{total_mass:.3f}-{mass_ratio:.3f}.dat"
            ),
            skip_header=1)
        return data[:,:2]

    except FileNotFoundError:
        raise ValueError({"error": "Requested data not found"})
