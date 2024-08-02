from typing import List
from scipy.optimize import curve_fit
import numpy as np


"""
    Transient Model Fitter

    The Transient Model Fitter defines models to estimate photometric
    magnitudes of observed astronomical objects. It is intended to be
    used in conjunction with the TypeScript based UI for the Transient
    Plotting Pool.

    There are two temporal models and two spectral models:
        - Temporal Models: Power Law and Exponential
        - Spectral Models: Power Law and Extinguished Power Law

    The `fit` method uses the `scipy.curve_fit` method to determine the
    best fitting parameters (temporal index, spectral index, and dust 
    extinction). See the `fit` method docstring for more details.
"""


class Model:
    def __init__(self, filters: List[str], ref_f: str, ref_t: int, ref_m: int, temporal: str):
        """ Photometric Model. Reference values are used to `anchor` the
        line to a starting point.

        :param filters: list of photometric filters
        :param ref_f: reference filter
        :param ref_t: reference time in units of days since trigger
        :param ref_m: reference magnitude
        :param temporal: temporal model, one of 'power law' or 'exponential'
        """
        self._FILTER_ZERO_POINT = {
            # Filter zero point values are empirically determined
            # See: https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html

            # Vega flux zero points in mJy (Bessell et al. (1998))
            'U': 1.79000,   'B': 4.06300,   'V': 3.63600,   'R': 3.06400,
            'I': 2.41600,   'J': 1.58900,   'H': 1.02100,   'K': 0.64000,

            # AB flux zero points in mJy (Fukugita et al. (1996))
            "u\'": 3.680,   "g\'": 3.643,   "r\'": 3.648,   "i\'": 3.644,
            "z\'": 3.631
        }

        self._FILTER_WAVELENGTH = {
            # Filter wavelengths in micro-meters
            'U': 0.36400,   'B': 0.44200,   'V': 0.54000,   'R': 0.64700,
            "I": 0.78650,   "u\'": 0.354,   "g\'": 0.475,   "r\'": 0.622,
            "i\'": 0.763,   "z\'": 0.905,   'J': 1.25000,   'H': 1.65000,
            'K': 2.15000,   'Ks': 2.1500,   'W1': 3.4000,   'W2': 4.6000,
            'W3': 12.000,   'W4': 22.000,   'BP': 0.5320,   'G': 0.67300,
            'RP': 0.7970
        }

        if temporal.lower() not in ['power law', 'exponential']:
            raise ValueError(f'Unsupported temporal model provided: {temporal}')

        self.data_filters = filters
        self.temporal = temporal.lower()
        self.ref_zero_point = self._FILTER_ZERO_POINT[ref_f]
        self.ref_wavelength = self._FILTER_WAVELENGTH[ref_f]
        self.ref_mag = ref_m
        self.ref_time = ref_t

    @staticmethod
    def calculate_extinction(filter_lambda: float) -> float:
        """ Calculates the average R_v dependent extinction factor.

        See: https://articles.adsabs.harvard.edu//full/1989ApJ...345..245C/0000249.000.html

        :param filter_lambda: the wavelength of the filter in meters
        :return: R_v dependent extinction factor
        """
        x = 1.0 / filter_lambda
        y = x - 1.82
        a, b = 0.0, 0.0

        # Infrared
        if 0.3 <= x <= 1.1:
            a = 0.574 * x ** 1.61
            b = -0.527 * x ** 1.61

        # Optical/NIR
        elif 1.1 < x <= 3.3:
            a = (1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 +
                 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.7753 * y ** 6 + 0.32999 * y ** 7)
            b = (1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 -
                 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.3026 * y ** 6 - 2.09002 * y ** 7)

        return a + b / 3.1

    def zero_points(self) -> List[float]:
        """ Returns a list of zero points corresponding to each filter.
        """
        return [self._FILTER_ZERO_POINT[f] for f in self.data_filters]

    def wavelengths(self) -> List[float]:
        """ Returns a list of wavelengths corresponding to each filter.
        """
        return [self._FILTER_WAVELENGTH[f] for f in self.data_filters]

    def get_zero_point_dependence(self) -> np.ndarray:
        """ Returns the zero point dependent factors for the magnitude
        calculation.
        """
        return np.log10(np.array(self.zero_points()) / self.ref_zero_point)

    def get_temporal_dependence(self, times: np.ndarray, index: float) -> np.ndarray:
        """ Calculates the temporal dependence for the model.

        :param times: array-like containing days since events
        :param index: temporal index
        :return: numpy array of temporal dependent factors
        """
        if self.temporal == 'power law':
            return np.log10((times / self.ref_time) ** index)

        return np.log10(np.exp(1)) * index * times / self.ref_time

    def get_spectral_dependence(self, index: float) -> np.ndarray:
        """ Calculates the spectral dependence for the model.

        :param index: spectral index
        :return: numpy array of spectral dependent factors
        """
        return np.log10((self.ref_wavelength / np.array(self.wavelengths())) ** index)

    def get_extinction_dependence(self, ebv: float) -> np.ndarray:
        """ Calculates the dust extinction dependence for the extinguished
        power law model.

        :param ebv: selective extinction
        :return: numpy array of extinction factors with A_v = 3.1.
        """
        return np.array([3.1 * self.calculate_extinction(wl) * ebv for wl in self.wavelengths()])

    def get_model(self, times: np.ndarray, a: float, b: float) -> np.ndarray:
        """ Models the photometric magnitudes.

        :param times: times in days since event
        :param a: temporal index
        :param b: spectral index
        :return: numpy array of modeled photometric magnitudes
        """
        eq_zp = self.get_zero_point_dependence()
        eq_time = self.get_temporal_dependence(times, a)
        eq_freq = self.get_spectral_dependence(b)

        return self.ref_mag - 2.5 * (eq_time + eq_freq - eq_zp)

    def get_extinction_model(self, times: np.ndarray, a: float, b: float, ebv: float):
        """ Models the photometric magnitudes.

        :param times: time in days since event
        :param a: temporal index
        :param b: spectral index
        :param ebv: e(b - v) extinction
        :return: numpy array of modeled photometric magnitudes
        """
        return self.get_model(times, a, b) + self.get_extinction_dependence(ebv)


def fit(xdata: List[float], ydata: List[float], filters: List[str], params) -> List[float]:
    """ Performs a non-linear regression fit on the provided data using
    scipy's curve_fit method. Returns the best fit parameters corresponding
    to the provided spectral model. E.g.,
        - power law model: temporal index, spectral index
        - extinguished power law model: temporal index, spectral index, E(B -V)

    To increase performance and fitting accuracy, initial guesses and bounds
    are used. Initial guesses are provided via the Transient Plotting Tool
    interface via the `params` request. The bounds are [generously] selected
    here based on typical values found for Gamma-Ray Burst (GRB) events. I.e.,
        - temporal index bound: [-3, 3] (negative indicates fading)
        - spectral index bound: [-3, 3]
        - ebv bounds: [0, 1]

    :param xdata: list of times in units of days since trigger
    :param ydata: list of photometric magnitudes
    :param filters: list of filters corresponding to the xdata and ydata
    :param params:
        - ref_f: reference filter
        - ref_t: reference time in units of days since trigger
        - ref_m: reference magnitude
        - ebv: dust extinction E(B-V) [used as initial guess]
        - a_index: temporal index (typically denoted by alpha) [used as initial guess]
        - b_index: spectral index (typically denoted by beta) [used as initial guess]
        - temporal: case-insensitive temporal model, one of 'power law' or 'exponential'
        - spectral: case-insensitive spectral model, one of 'power law' or 'extinguished power law'
    :return: list of best fit parameters rounded to 3 sig figs
    """

    # Use numpy arrays for operational (multiplication, division) convenience
    xdata, ydata = np.array(xdata), np.array(ydata)

    # Initialize the photometric magnitude model
    model = Model(filters, params['ref_f'], params['ref_t'], params['ref_m'], params['temporal'])

    # Model the data to determine the best temporal and spectral indices
    if params['spectral'].lower() == 'power law':
        popt = curve_fit(model.get_model, xdata, ydata,
                         p0=[params['a_index'], params['b_index']],
                         bounds=([-3, -3], [3, 3]))

    # Model the data to determine the best temporal index, spectral index, and E(B-V)
    elif params['spectral'].lower() == 'extinguished power law':
        popt = curve_fit(model.get_extinction_model, xdata, ydata,
                         p0=[params['a_index'], params['b_index'], params['ebv']],
                         bounds=([-3, -3, 0], [3, 3, 1]))

    else:  # Invalid model provided
        raise ValueError(f"Unsupported spectral model provided: {params['spectral']}")

    return [round(p, 3) for p in popt[0]]  # Best fit values
