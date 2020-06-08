"""

ariel@padova
02/06/2020

Miscelaneous tools.

"""

import numpy as np


def calc_sn(wl, f_obs, f_err, z, window_limits=(5500, 5700)):

    wl = wl / (1+z)

    window = (wl > window_limits[0]) & (wl < window_limits[1])

    sn = np.mean(f_obs[window]) / np.mean(np.sqrt(f_err[window]))

    return sn


def calc_continuum_rms(wl, f_obs, f_syn, z, window_limits=(5500, 5700)):

    wl = wl / (1+z)

    window = (wl > window_limits[0]) & (wl < window_limits[1])

    continuum_rms = np.mean(np.absolute(f_obs[window] - f_syn[window]))

    return continuum_rms

