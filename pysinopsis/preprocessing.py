import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage.filters import gaussian
from skimage.transform import rescale
from scipy.interpolate import interp1d
import itertools
from pysinopsis.utils import resampler


def smooth_cube(cube, sigma, gauss_truncate=5):
    """
    Smooth a cube using a Gaussian filter.

    Parameters:
        cube (numpy.ndarray): Input cube to be smoothed.
        sigma (float): Standard deviation for Gaussian kernel. 
        gauss_truncate (float, optional): Truncate parameter for Gaussian kernel. Default is 5.

    Returns:
        numpy.ndarray: Smoothed cube.
    """

    smoothed_cube = np.empty_like(cube)

    for i in range(cube.shape[0]):
        smoothed_cube[i, :, :] = gaussian(cube[i, :, :], sigma, truncate=gauss_truncate)

    return smoothed_cube


def bin_cube(cube, bin_size):
    """
    Bin a cube by reducing its spatial dimensions.

    Parameters:
        cube (numpy.ndarray): Input cube to be binned.
        bin_size (int): Size of the binning factor. The cube will be reduced by a factor of `bin_size`
            along each spatial dimension.

    Returns:
        numpy.ndarray: Binned cube.
    """

    new_scale = 1 / bin_size

    binned_cube = np.empty(shape=(cube.shape[0], int(cube.shape[1]/bin_size), int(cube.shape[2]/bin_size)))

    for i in range(cube.shape[0]):
        binned_cube[i, :, :] = rescale(cube[i, :, :].astype(float), new_scale, anti_aliasing=False)

    return binned_cube


def process_degraded_cube(in_filename, out_filename, psf_sigma, bin_size, gauss_truncate=5):
    """
    Degrades cube by smoothing and binning.

    Parameters:
        in_filename (str): Filename of the input cube FITS file.
        out_filename (str): Filename for the output processed cube FITS file.
        psf_sigma (float): Standard deviation of the Gaussian PSF used for smoothing.
        bin_size (int): Size of the binning factor for spatial dimensions.
        gauss_truncate (int, optional): Truncation factor for the Gaussian kernel used in smoothing. Defaults to 5.

    Returns:
        HDUList: FITS HDUList object containing the processed cube.
    """

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


def integrated_spectrum(wl, flux_cube, error_cube, z, mask, wl_range):
    """
    Computes the integrated spectrum over a given wavelength range.

    Parameters:
        wl (array-like): Wavelength array.
        flux_cube (array-like): 3D array of flux values.
        error_cube (array-like): 3D array of error values.
        z (array-like): Redshift array.
        mask (array-like): 3D boolean mask array.
        wl_range (tuple): Tuple representing the range of wavelengths to integrate over.

    Returns:
        tuple: Tuple containing the resampled wavelength array, integrated flux array, and integrated error array.
    """

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
