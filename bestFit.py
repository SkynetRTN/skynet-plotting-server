
try:
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    from time import time
    import numpy as np
except ImportError as e:
    print(e)


class InitialGuess:
    def __init__(self, m, a, b):
        self.m = m
        self.a = a
        self.b = b


class Model:
    def __init__(self, ref_fltr, fltrs, x):
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

        self.data_fltrs = fltrs
        self.ref_zero_point = self._FILTER_ZERO_POINT[ref_fltr]
        self.ref_wavelength = self._FILTER_WAVELENGTH[ref_fltr]
        self.referenceX = x

    def _zero_point_list(self, d):
        return [self._FILTER_ZERO_POINT[d[key]] for key in d.keys()]

    def _wavelength_list(self, d):
        return [self._FILTER_WAVELENGTH[d[key]] for key in d.keys()]       

    def get_model(self, x, m, a, b):
        zero_point = np.array(self._zero_point_list(self.data_fltrs))
        wavelength = np.array(self._wavelength_list(self.data_fltrs))

        eq_zp = np.log10(self.ref_zero_point / zero_point)
        eq_time = a * np.log10(x / self.referenceX)
        eq_freq = b * np.log10(wavelength / self.ref_wavelength)

        return m - 2.5 * (eq_zp + eq_time + eq_freq)


def getBestFit(model, x, y, guess):
    return curve_fit(model, x, y, \
        p0=[guess.m, guess.a, guess.b], \
        bounds=([-10, -3, -2], [20, 1, 1]))


def fitToData(xdata, ydata, filters, params):

    guess = InitialGuess(params['m'], params['a'], params['b'])
    model = Model(params['filter'], filters, params['t'])

    # optimal, covariance
    popt, pcov = getBestFit(model.get_model, xdata, ydata, guess)
    return [round(p, 3) for p in popt]