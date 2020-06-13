"""

ariel@padova
29/05/2020

Tools to organize SINOPSIS output in python.

TODO: Method for radial profiles, requires calculating distance to centre of the cube (wcs?)
TODO: Past and future datacubes, require parametric SFH module and spectral library handling

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from numpy.ma import masked_array
from astropy.io import fits
from astropy.table import Table
from collections import OrderedDict
import pysinopsis.plotting as sinplot


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

    def __init__(self, galaxy_id, obs_file, sinopsis_directory='./'):

        self.galaxy_id = galaxy_id
        self.obs_file = obs_file
        self.sinopsis_directory = sinopsis_directory

        # Ages:
        self.age_bins = np.genfromtxt(sinopsis_directory + galaxy_id + '.log', skip_header=22, skip_footer=6)[:, 0]

        # SINOPSIS results:
        self.header_info, self.properties = read_results_cube(sinopsis_directory + galaxy_id + '_out.fits')

        # Equivalend widths:
        self.eqw = read_eqw_cube(sinopsis_directory + galaxy_id + '_eqw.fits')

        # Magnitudes:
        self.mag = read_mag_cube(sinopsis_directory + galaxy_id + '_mag.fits')

        # Observed cube:
        obs_cube = fits.open(obs_file)

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

    def invalid_spaxel(self, x, y):
        if np.all(self.f_syn.mask[:, x, y]):
            return True
        else:
            return False

    def fit_details(self, x, y):

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        else:

            fname_prefix = self.obs_file.split('.')[0]
            fname_fit_details = fname_prefix + '.%0.4d_%0.4d.out' % (x+1, y+1)  # FIXME: 0- or 1-indexed?

            n_features = open(fname_fit_details, 'rt').read().splitlines()[0].split()[0]
            n_features = int(n_features)

            fit_details_table = Table.read(fname_fit_details, format='ascii', header_start=1,
                                           data_end=n_features+2)

            return fit_details_table

    def plot_map(self, sinopsis_property):

        plt.figure()

        sinplot.plot_sinopsis_map(self, sinopsis_property)

        plt.show()

    def plot_spectrum(self, x, y, plot_error=True, plot_legend=True):

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        plt.figure()

        # FIXME : BUNIT conversion sort of hard-coded
        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y], self.f_err[:, x, y],
                         plot_error=plot_error, plot_legend=plot_legend,
                         flux_unit=float(self.obs_header['primary']['BUNIT_JA'][0: 5]))

        plt.show()

    def plot_fit_complete(self, x, y, figsize=(7.75, 6.5)):

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        plt.figure(figsize=figsize)

        gs1 = gridspec.GridSpec(3, 3)
        gs1.update(bottom=0.47, top=0.96, hspace=0.05, right=0.96)

        gs2 = gridspec.GridSpec(1, 3)
        gs2.update(top=0.38, bottom=0.08, hspace=0.05, right=0.96, wspace=0.02)

        ax_spectrum = plt.subplot(gs1[0:2, :])
        ax_residuals = plt.subplot(gs1[2, :], sharex=ax_spectrum)
        ax_sfh = plt.subplot(gs2[0, 0])
        ax_annotations = plt.subplot(gs2[0, 1:3])

        # FIXME : BUNIT conversion sort of hard-coded
        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y], self.f_err[:, x, y],
                         flux_unit=float(self.obs_header['primary']['BUNIT_JA'][0: 5]), ax=ax_spectrum)

        sinplot.plot_residuals(self.wl, self.f_obs[:, x, y], self.f_syn_cont[:, x, y], ax=ax_residuals)


if __name__ == '__main__':
    import importlib
    importlib.reload(sinplot)

    sinopsis_cube = SinopsisCube('A370_01', 'tests/test_run/A370_01_DATACUBE_FINAL_v1_ec.fits', 'tests/test_run/')
    # sinopsis_cube.plot_spectrum(25, 27, plot_error=False)
    # sinopsis_cube.plot_map('SFR1')
    sinopsis_cube.plot_fit_complete(25, 27)