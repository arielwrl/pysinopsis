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
from pysinopsis import utils


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
    for i in range(77):
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
        self.age_bins = np.append(0, self.age_bins)
        self.age_bins_4 = np.genfromtxt(sinopsis_directory + galaxy_id + '.bin', skip_header=3)
        self.age_bins_4 = np.append(0, self.age_bins_4)

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

    def plot_map(self, sinopsis_property, show_plot=True):

        plt.figure()

        sinplot.plot_sinopsis_map(self, sinopsis_property)

        if show_plot:
            plt.show()

    def plot_spectrum(self, x, y, plot_error=True, plot_legend=True, show_plot=True):

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        plt.figure()

        # FIXME : BUNIT conversion sort of hard-coded
        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y], self.f_err[:, x, y],
                         plot_error=plot_error, plot_legend=plot_legend,
                         flux_unit=float(self.obs_header['primary']['BUNIT_JA'][0: 5]))

        if show_plot:
            plt.show()

    def plot_fit_complete(self, x, y, figsize=(13.25, 8.5), out_file_name=None, out_format='png',
                          show_plot=True):

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        fig = plt.figure(figsize=figsize)

        gs = gridspec.GridSpec(5, 5)

        ax_spectrum = plt.subplot(gs[0:2, 0:5])
        ax_residuals = plt.subplot(gs[2, 0:5], sharex=ax_spectrum)
        ax_sfh = plt.subplot(gs[3:5, 0:2])
        ax_map = plt.subplot(gs[3:5, 2:4])

        # Plotting fit and residuals:
        # FIXME: BUNIT conversion sort of hard-coded
        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y], self.f_err[:, x, y],
                         flux_unit=float(self.obs_header['primary']['BUNIT_JA'][0: 5]), ax=ax_spectrum)

        sinplot.plot_residuals(self.wl, self.f_obs[:, x, y], self.f_syn_cont[:, x, y], ax=ax_residuals)

        # Plotting SFH:
        sfr = np.array([self.properties['sfr_'+str(i)][x, y] for i in range(1, 13)])
        sinplot.plot_sfh(self.age_bins, sfr, ax=ax_sfh)

        # A map to show the spaxel:
        # FIXME: Plotting map of total flux, did not think enough about this
        # FIXME: Inverting x and y !!! Have to find a better solution for this
        reference_map = ax_map.imshow(np.sum(self.f_obs, axis=0), cmap='Blues', origin='lower')
        fig.colorbar(mappable=reference_map, ax=ax_map, label=r'$\log F_\lambda$')
        ax_map.scatter(y, x, marker='x', s=80, c='r', label='x =' + str(y) + ', y =' + str(x))
        ax_map.legend()
        ax_map.set_xlabel('x', fontsize=12)
        ax_map.set_ylabel('y', fontsize=12)

        # Annotating results:
        # FIXME: This code is hideous
        property_list = ['lwage', 'mwage', 'Av', 'Av_y']
        property_values = [self.properties[sinopsis_property][x, y] for sinopsis_property in property_list]
        uncertainty_list = [utils.get_uncertainty(self, sinopsis_property, x, y) for sinopsis_property in property_list]
        annotations = [sinplot.sinopsis_labels[property_list[i]] + r'$ = %0.2f^{+%0.2f}_{-%0.2f}$'
                       % (property_values[i], uncertainty_list[i][0], uncertainty_list[i][1])
                       for i in range(len(property_list))]

        for i in range(len(annotations)):
            ax_residuals.annotate(annotations[i], xy=(0.82, 0.3-i*0.05), xycoords='figure fraction', fontsize=12)

        fig.subplots_adjust(top=0.98,
                            bottom=0.06,
                            left=0.06,
                            right=0.99,
                            hspace=0.389,
                            wspace=0.213)

        if out_file_name is not None:
            plt.savefig(out_file_name, format=out_format)

        if show_plot:
            plt.show()
        else:
            plt.close()


if __name__ == '__main__':
    import importlib
    importlib.reload(sinplot)
    importlib.reload(utils)

    sinopsis_cube = SinopsisCube('A2744_06', 'tests/test_run/A2744_06_DATACUBE_FINAL_v1_ec.fits', 'tests/test_run/')
    # sinopsis_cube.plot_spectrum(25, 27, plot_error=False)
    # sinopsis_cube.plot_map('SFR1')
    sinopsis_cube.plot_fit_complete(54, 29)

