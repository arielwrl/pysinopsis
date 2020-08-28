"""

ariel@padova
28/08/2020

Tools to perform several pre-processing steps

"""

import numpy as np
import matplotlib.pyplot as plt
import itertools
from pysinopsis.output import SinopsisCube
from skimage.filters import gaussian


def cube_smooth(cube, sigma, gauss_truncate=5):
    """

    :param cube: masked_array
    flux, i, j
    :param sigma:
    :param gauss_truncate:
    :return:

    """

    smoothed_cube = np.empty_like(cube)

    for i in range(cube.shape[0]):
        smoothed_cube[i, :, :] = gaussian(cube[i, :, :], sigma, truncate=gauss_truncate)

    return smoothed_cube


if __name__ == '__main__':

    sinopsis_cube = SinopsisCube('A2744_06', 'tests/test_run/A2744_06_DATACUBE_FINAL_v1_ec.fits', 'tests/test_run/')

    processed_cube = cube_smooth(sinopsis_cube.f_obs, 5)

    plt.figure()
    plt.imshow(processed_cube[1000, :, :], origin='lower')

    plt.figure()
    plt.plot(sinopsis_cube.wl, sinopsis_cube.f_obs[:, 54, 29])
    plt.plot(sinopsis_cube.wl, processed_cube[:, 54, 29])

