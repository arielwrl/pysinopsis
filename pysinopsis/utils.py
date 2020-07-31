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


def get_uncertainty(sinopsis_cube, sinopsis_property, x, y):

    if sinopsis_property in ['Mb1', 'Mb2', 'Mb3', 'Mb4']:
        print('Uncertainty in masses to be automated.')

    else:
        plus = sinopsis_cube.properties[sinopsis_property+'_M'][x, y] - \
               sinopsis_cube.properties[sinopsis_property][x, y]
        minus = sinopsis_cube.properties[sinopsis_property][x, y] - \
                sinopsis_cube.properties[sinopsis_property+'_m'][x, y]

        return plus, minus

