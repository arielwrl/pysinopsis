"""

ariel@padova
02/06/2020

Miscelaneous tools.

"""

import numpy as np


def calc_sn(wl, f_obs, f_err, z, window_limits=(5500, 5700)):
    wl = wl / (1 + z)

    window = (wl > window_limits[0]) & (wl < window_limits[1])

    sn = np.mean(f_obs[window]) / np.mean(np.sqrt(f_err[window]))

    return sn


def calc_continuum_rms(wl, f_obs, f_syn, z, window_limits=(5500, 5700)):
    wl = wl / (1 + z)

    window = (wl > window_limits[0]) & (wl < window_limits[1])

    continuum_rms = np.mean(np.absolute(f_obs[window] - f_syn[window]))

    return continuum_rms


def get_uncertainty(sinopsis_cube, sinopsis_property, x, y):
    if sinopsis_property in ['Mb1', 'Mb2', 'Mb3', 'Mb4']:
        print('Uncertainty in masses to be automated.')

    else:
        plus = sinopsis_cube.properties[sinopsis_property + '_M'][x, y] - \
               sinopsis_cube.properties[sinopsis_property][x, y]
        minus = sinopsis_cube.properties[sinopsis_property][x, y] - \
                sinopsis_cube.properties[sinopsis_property + '_m'][x, y]

        return plus, minus


def initial_burst(t, t_u, n1, tau_i):
    """

    Equation (2) in page 29 of SINOPSIS manual, version 1.6.7

    """

    sfr = (((t_u - t) / t_u) ** n1) * np.exp(-((t_u - t) / (tau_i * t_u)))

    return sfr


def late_burst(t, m_b, t_b, n2, tau_b):
    sfr = m_b * (((t_b - t) / t_b) ** n2) * np.exp(-((t_b - t) / (tau_b * t_b)))

    return sfr


def get_center(fit_mask):
    coords = np.argwhere(fit_mask == 1)

    x_center = np.mean(coords[:, 1])
    y_center = np.mean(coords[:, 0])

    return x_center, y_center


def get_coordinate_crid(sinopsis_cube):
    x_coords, y_coords = np.meshgrid(range(sinopsis_cube.cube_shape[1]), range(sinopsis_cube.cube_shape[2]),
                                     indexing='xy')
    grid = np.array([x_coords, y_coords])

    return grid


def gini(x):
    # Mean absolute difference
    mad = np.abs(np.subtract.outer(x, x)).mean()
    # Relative mean absolute difference
    rmad = mad/np.mean(x)
    # Gini coefficient
    g = 0.5 * rmad
    return g



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    t = np.linspace(1e6, 1.4e9, 1000)
    initial = initial_burst(t, t_u=1.4e9, n1=0.01, tau_i=0.05)
    plt.plot(np.log10(t), initial)

    t = np.linspace(1.e6, 1.e8, 1000)
    burst = late_burst(t, m_b=0.1, t_b=1.0e8, n2=0.1, tau_b=1)
    plt.plot(np.log10(t), burst)
