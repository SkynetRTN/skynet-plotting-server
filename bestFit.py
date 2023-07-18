
from typing import List


try:
    from scipy.optimize import curve_fit
    import numpy as np
except ImportError as e:
    print(e)


class InitialGuess:
    def __init__(self, m, a, b):
        self.m = m
        self.a = a
        self.b = b


class Model:
    def __init__(self, ref_filter, filters, x, m):
        self._FILTER_ZERO_POINT = {
            'U' : 1.790,
            'B' : 4.063,
            'V' : 3.636,
            'R' : 3.064,
            'I' : 2.416,
            'J' : 1.589,
            'H' : 1.021,
            'K' : 0.640,
            "u\'": 3.680,
            "g\'": 3.643,
            "r\'": 3.648,
            "i\'": 3.644,
            "z\'": 3.631,
        }

        self._FILTER_WAVELENGTH = {
            'U': 0.364,
            'B': 0.442,
            'V': 0.54,
            'R': 0.647,
            "I": 0.7865,
            "u\'": 0.354,
            "g\'": 0.475,
            "r\'": 0.622,
            "i\'": 0.763,
            "z\'": 0.905,
            'J': 1.25,
            'H': 1.65,
            'K': 2.15,
            'Ks': 2.15,
        }

        self.data_filters = filters
        self.ref_zero_point = self._FILTER_ZERO_POINT[ref_filter]
        self.ref_wavelength = self._FILTER_WAVELENGTH[ref_filter]
        self.referenceX = x
        self.referenceMag = m

    def zero_point_list(self, filters: List[str]):
        return [self._FILTER_ZERO_POINT[f] for f in filters]
        # return [self._FILTER_ZERO_POINT[d[key]] for key in d.keys()]

    def wavelength_list(self, filters: List[str]):
        return [self._FILTER_WAVELENGTH[f] for f in filters]
        # return [self._FILTER_WAVELENGTH[d[key]] for key in d.keys()]

    def get_model(self, x, a, b):
        zero_point = np.array(self.zero_point_list(self.data_filters))
        wavelength = np.array(self.wavelength_list(self.data_filters))

        eq_zp = np.log10(self.ref_zero_point / zero_point)
        eq_time = a * np.log10(x / self.referenceX)
        eq_freq = b * np.log10(wavelength / self.ref_wavelength)
        # extinction term

        return self.referenceMag - 2.5 * (eq_zp + eq_time + eq_freq)


def getBestFit(model, x, y, guess):
    return curve_fit(model, x, y, \
        p0=[guess.a, guess.b], \
        bounds=([-3, -2], [1, 1]))


def fitToData(xdata: List[float], ydata: List[float], filters, params):

    try:
        guess = InitialGuess(params['m'], params['a'], params['b'])
        model = Model(params['filter'], filters, params['t'], params['m'])
    except Exception as e:
        print(e)
        return []

    # optimal, covariance
    popt, pcov = getBestFit(model.get_model, xdata, ydata, guess)
    return [round(p, 3) for p in popt]