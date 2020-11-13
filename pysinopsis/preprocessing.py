"""

ariel@padova
28/08/2020

Tools to perform several pre-processing steps

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage.filters import gaussian
from skimage.transform import rescale
from scipy.interpolate import interp1d
import itertools


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


def process_degraded_cube(in_filename, out_filename, psf_sigma, bin_size, gauss_truncate=5):

    cube = fits.open(in_filename)

    smoothed_flux = smooth_cube(cube[1].data, psf_sigma, gauss_truncate)
    smoothed_error = smooth_cube(cube[2].data, psf_sigma, gauss_truncate)

    binned_flux = bin_cube(smoothed_flux, bin_size)
    binned_error = bin_cube(smoothed_error, bin_size)

    hdu_list = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=binned_flux), fits.ImageHDU(data=binned_error)])

    for i in range(3):
        hdu_list[i].header = cube[i].header

    hdu_list.writeto(out_filename)

    return hdu_list


def resampler(x_old, y_old, x_new):
    interp = interp1d(x_old, y_old, bounds_error=False, fill_value=(0., 0.))

    y_new = interp(x_new)

    return y_new


def integrated_spectrum(wl, flux_cube, error_cube, z, mask, wl_range):

    flux_cube = np.ma.masked_array(flux_cube, mask=mask)
    error_cube = np.ma.masked_array(error_cube, mask=mask)

    z = z - np.mean(z)

    wl_resampled = np.arange(wl_range[0], wl_range[1], 3).astype(float)

    integrated_flux = np.zeros_like(wl_resampled)
    integrated_error = np.zeros_like(wl_resampled)

    for i, j in itertools.product(range(flux_cube.data.shape[1]), range(flux_cube.data.shape[2])):

        if np.any(mask[:, i, j]):
            continue

        wl_rest = wl / (1 + z[i, j])

        flux_interp = resampler(wl_rest, flux_cube[:, i, j], wl_resampled)
        error_interp = resampler(wl_rest, error_cube[:, i, j], wl_resampled)

        integrated_flux += flux_interp
        integrated_error += error_interp

    return wl_resampled, integrated_flux, integrated_error


if __name__ == '__main__':

    degraded_cube = process_degraded_cube('tests/test_run/A2744_06_DATACUBE_FINAL_v1_ec.fits',
                                          'tests/A2744_06_DATACUBE_FINAL_v1_ec_DEGRADED.fits', psf_sigma=5, bin_size=2)

    # sinopsis_cube = SinopsisCube('tests/test_run/')
    #
    # smoothed = smooth_cube(sinopsis_cube.f_obs, 5)
    #
    # plt.figure()
    # plt.imshow(smoothed[1000, :, :], origin='lower')
    #
    # plt.figure()
    # plt.plot(sinopsis_cube.wl, sinopsis_cube.f_obs[:, 54, 29])
    # plt.plot(sinopsis_cube.wl, smoothed[:, 54, 29])
    #
    # binned = bin_cube(smoothed, 10)
    #
    # plt.figure()
    # plt.imshow(binned[1000, :, :], origin='lower')
    #
    # plt.figure()
    # plt.plot(sinopsis_cube.wl, sinopsis_cube.f_obs[:, 54, 29])
    # plt.plot(sinopsis_cube.wl, binned[:, 6, 2])
