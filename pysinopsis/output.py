"""

ariel@padova
29/05/2020

Tools to organize SINOPSIS output in python.

TODO: Method for radial profiles, requires calculating distance to centre of the cube (wcs?)
TODO: Past and future datacubes, require parametric SFH module and spectral library handling
TODO: Difference between observed and synthetic datacube

"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.ma import masked_array
from astropy.io import fits
from collections import OrderedDict
from pysinopsis.plotting import plot_fit, plot_sinopsis_map


def read_results_cube(output_cube_file):
    """

    Reads SINOPSIS output cube into a dictionary of masked arrays.

    input:
    -----------
    output_cube_file: file name of SINOPSIS output cube
        type: str

    returns:
    -----------
    dictionary containing SINOPSIS outputs

    """

    output_cube = fits.open(output_cube_file)[0]

    header_info_keys = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'EXTEND', 'CRVAL1', 'CRVAL2',
                        'CD1_1', 'CD2_1', 'CD1_2', 'CD2_2', 'CRPIX1', 'CRPIX2', 'CTYPE1', 'CTYPE2', 'OBJECT']

    header_info = OrderedDict()
    for key in header_info_keys:
        header_info[key] = output_cube.header[key]

    properties = OrderedDict()
    for i in range(76):
        plane_key = 'PLANE%0.2d' % i

        properties[output_cube.header[plane_key]] = masked_array(output_cube.data[i], mask=output_cube.data[i] == -999)

    return header_info, properties


def read_eqw_cube(eqw_cube_file):
    """

    Reads SINOPSIS eqw cube into a dictionary of masked arrays.

    input:
    -----------
    output_cube_file: file name of SINOPSIS equivalent widths cube
        type: str

    returns:
    -----------
    dictionary containing SINOPSIS eqw

    """

    eqw_cube = fits.open(eqw_cube_file)[0]

    eqws = OrderedDict()
    for i in range(31):
        plane_key = 'PLANE%0.2d' % i  # FIXME: Hard-coded!

        eqws[eqw_cube.header[plane_key]] = masked_array(eqw_cube.data[i], mask=eqw_cube.data[i] == -999)

    return eqws


def read_mag_cube(mag_cube_file):
    """

    Reads SINOPSIS magnitudes cube into a dictionary of masked arrays.

    input:
    -----------
    output_cube_file: file name of SINOPSIS equivalent widths cube
        type: str

    returns:
    -----------
    dictionary containing SINOPSIS magnitudes

    """

    mag_cube = fits.open(mag_cube_file)[0]

    mags = OrderedDict()
    for i in range(33):
        plane_key = 'PLANE%0.2d' % i  # FIXME: Hard-coded!

        mags[mag_cube.header[plane_key]] = masked_array(mag_cube.data[i], mask=mag_cube.data[i] == -999)

    return mags


class SinopsisCube:
    """

    Reads SINOPSIS output for a datacube into a python object.

    input:
    -----------
    sinopsis_directory: Directory with SINOPSIS files (must end with /)
        type: str

    galaxy_id: ID of the galaxy, used to read output files
        type: str

    """

    def __init__(self, galaxy_id, obs_file, sinopsis_directory):

        # Ages:
        self.age_bins = np.genfromtxt(sinopsis_directory + galaxy_id + '.log', skip_header=22, skip_footer=6)[:, 0]

        # SINOPSIS results:
        self.header_info, self.properties = read_results_cube(sinopsis_directory + galaxy_id + '_out.fits')

        # Equivalend widths:
        self.eqw = read_eqw_cube(sinopsis_directory + galaxy_id + '_eqw.fits')

        # Magnitudes:
        self.mag = read_mag_cube(sinopsis_directory + galaxy_id + '_mag.fits')

        # Observed cube:
        obs_cube = fits.open(sinopsis_directory + obs_file)

        self.obs_header = {'primary': obs_cube[0].header,
                           'flux': obs_cube[1].header,
                           'error': obs_cube[2].header}

        self.wl = obs_cube[1].header['CRVAL3'] + obs_cube[1].header['CD3_3'] * np.arange(obs_cube[1].header['NAXIS3'])
        self.f_obs = masked_array(obs_cube[1].data, mask=np.isnan(obs_cube[1].data))
        self.f_err = masked_array(np.sqrt(obs_cube[2].data), mask=np.isnan(obs_cube[2].data))  # FIXME: CHECK!
        self.cube_shape = self.f_obs.shape

        # Model cube:
        model_cube = fits.open(sinopsis_directory + galaxy_id + '_modelcube.fits')[0]
        model_cube_nolines = fits.open(sinopsis_directory + galaxy_id + '_modelcube_nolines.fits')[0]

        self.f_syn = masked_array(model_cube.data, mask=model_cube.data == -999)
        self.f_syn_cont = masked_array(model_cube_nolines.data, mask=model_cube.data == -999)

        # Emission only cube:
        self.emission_only = self.f_obs - self.f_syn_cont
        self.emission_only_model = self.f_syn - self.f_syn_cont  # FIXME: Does it make sense to have this?

    def plot_spaxel(self, x, y, plot_error=True, plot_legend=True):

        plt.figure()

        # FIXME : BUNIT conversion sort of hard-coded
        plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y], self.f_err[:, x, y],
                 plot_error=plot_error, plot_legend=plot_legend,
                 flux_unit=float(self.obs_header['primary']['BUNIT_JA'][0: 5]))

        plt.show()

    def plot_map(self, sinopsis_property):

        plt.figure()

        plot_sinopsis_map(self, sinopsis_property)

        plt.show()


if __name__ == '__main__':
    sinopsis_cube = SinopsisCube('A370_01', 'A370_01_DATACUBE_FINAL_v1_ec.fits', 'tests/test_run/')
    sinopsis_cube.plot_spaxel(25, 25, plot_error=False)
    sinopsis_cube.plot_map('SFR1')

