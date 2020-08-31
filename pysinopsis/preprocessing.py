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
from skimage.transform import rescale


def smooth_cube(cube, sigma, gauss_truncate=5):
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


def bin_cube(cube, bin_size):
    """

    :param cube:
    :param bin_size:
        number of pixels that will be binned together
    :return:
    """

    new_scale = 1 / bin_size

    binned_cube = np.empty(shape=(cube.shape[0], int(cube.shape[1]/bin_size), int(cube.shape[2]/bin_size)))

    for i in range(cube.shape[0]):
        binned_cube[i, :, :] = rescale(cube[i, :, :].astype(float), new_scale, anti_aliasing=False)
    # Note: .astype(float) is a workaround for what seems to be a bug in scikit-image

    return binned_cube


if __name__ == '__main__':

    sinopsis_cube = SinopsisCube('A2744_06', 'tests/test_run/A2744_06_DATACUBE_FINAL_v1_ec.fits', 'tests/test_run/')

    smoothed = smooth_cube(sinopsis_cube.f_obs, 5)

    plt.figure()
    plt.imshow(smoothed[1000, :, :], origin='lower')

    plt.figure()
    plt.plot(sinopsis_cube.wl, sinopsis_cube.f_obs[:, 54, 29])
    plt.plot(sinopsis_cube.wl, smoothed[:, 54, 29])

    binned = bin_cube(smoothed, 10)

    plt.figure()
    plt.imshow(binned[1000, :, :], origin='lower')

    plt.figure()
    plt.plot(sinopsis_cube.wl, sinopsis_cube.f_obs[:, 54, 29])
    plt.plot(sinopsis_cube.wl, binned[:, 6, 2])
