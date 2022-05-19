import numpy as np


def find_gravity_data(mass_ratio, total_mass):
    try:
        data = np.genfromtxt(
            os.path.join(
                os.path.dirname(__file__),
                "gravity_dat_data", f"mt_{total_mass}"
                f"simulation-TD-SEOBNRv2-{total_mass:.5f}-{larger_mass:.5f}.dat",
            ),
            skip_header=1)

    except FileNotFoundError:
        raise ValueError({"error": "Requested data not found"})