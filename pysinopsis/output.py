"""

ariel@padova
29/05/2020

Tools to organize SINOPSIS output in python.

TODO: Method for radial profiles, requires calculating distance to centre of the cube (wcs?)

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


def read_config(sinopsis_dir='./'):
    """

    :param sinopsis_dir:
    :return: config dictionary
    """

    config_file = np.array(open(sinopsis_dir + 'config.sin', 'rb').readlines()).astype(str)

    config_dict = {'input_catalog': config_file[12].split()[-1],
                   'input_type': config_file[19].split()[-1],
                   'sinopsis_dir': sinopsis_dir,
                   'sfh_type': config_file[59].split()[-1],
                   }

    if config_dict['input_type'] == 'cube':
        config_dict['galaxy_id'] = config_dict['input_catalog'].split('.')[0]

    return config_dict


def read_sinopsis_catalog(sinopsis_dir, input_catalog_file, input_type):
    """

    :param sinopsis_dir:
    :param input_catalog_file:
    :param input_type:
    :return:
    """

    if input_type == 'cube':

        input_catalog = np.array(open(sinopsis_dir + input_catalog_file, 'rb').readlines()).astype(str)

        sinopsis_catalog = {'obs_file': input_catalog[0].split()[0],
                            'abs_mask_file': input_catalog[1].split()[0],
                            'z': float(input_catalog[2])
                            }
        try:
            sinopsis_catalog['em_mask_file'] = input_catalog[1].split()[1]
        except Exception:
            print('No emission line mask')

        return sinopsis_catalog


def read_results_cube(output_cube_file, sfh_type):
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

    global len_cube
    output_cube = fits.open(output_cube_file)[0]

    header_info_keys = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'EXTEND', 'CRVAL1', 'CRVAL2',
                        'CD1_1', 'CD2_1', 'CD1_2', 'CD2_2', 'CRPIX1', 'CRPIX2', 'CTYPE1', 'CTYPE2', 'OBJECT']

    header_info = OrderedDict()
    for key in header_info_keys:
        header_info[key] = output_cube.header[key]

    properties = OrderedDict()

    if sfh_type == 'ff':
        len_cube = 89
    if sfh_type == 'dexp':
        len_cube = 74

    for i in range(len_cube):
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
    for i in range(eqw_cube.data.shape[0]):
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
    for i in range(mag_cube.data.shape[0]):
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

    """

    def __init__(self, sinopsis_directory='./'):

        config = read_config(sinopsis_directory)
        catalog = read_sinopsis_catalog(sinopsis_directory, config['input_catalog'], input_type='cube')

        self.config = config
        self.catalog = catalog

        self.galaxy_id = config['galaxy_id']
        self.obs_file = catalog['obs_file']
        self.sinopsis_directory = sinopsis_directory

        # SINOPSIS results:
        self.header_info, self.properties = read_results_cube(sinopsis_directory + self.galaxy_id + '_out.fits',
                                                              sfh_type=self.config['sfh_type'])

        # Ages:
        if self.config['sfh_type'] == 'ff':
            self.age_bins = np.genfromtxt(sinopsis_directory + self.galaxy_id + '.log', skip_header=22,
                                          skip_footer=6)[:, 0]  # FIXME: Double-check me
            self.age_bins = np.append(0, self.age_bins)
            self.age_bins_4 = np.genfromtxt(sinopsis_directory + self.galaxy_id + '.bin', skip_header=3)
            self.age_bins_4 = np.append(0, self.age_bins_4)

            self.age_bin_center = np.array([(self.age_bins[i] + self.age_bins[i + 1]) / 2 for i in range(12)])

            self.sfh = np.array([self.properties['sfr_' + str(i)] for i in range(1, 13)])
            self.sfh = np.ma.masked_array(self.sfh, mask=self.sfh == -999)

        # Equivalend widths:
        self.eqw = read_eqw_cube(sinopsis_directory + self.galaxy_id + '_eqw.fits')

        # Magnitudes:
        self.mag = read_mag_cube(sinopsis_directory + self.galaxy_id + '_mag.fits')

        # Mask:
        self.mask = fits.open(sinopsis_directory + self.galaxy_id + '_fitmask.fits')[0].data

        # Observed cube:
        obs_cube = fits.open(self.sinopsis_directory + self.obs_file)

        self.obs_header = {'primary': obs_cube[0].header,
                           'flux': obs_cube[1].header,
                           'error': obs_cube[2].header}

        self.flux_unit = 10 ** float(self.obs_header['flux']['BUNIT'].split()[0].split('**')[1])
        self.flux_unit_err = 10 ** float(self.obs_header['error']['BUNIT'].split()[0].split('**')[1])

        self.wl = obs_cube[1].header['CRVAL3'] + obs_cube[1].header['CD3_3'] * np.arange(obs_cube[1].header['NAXIS3'])
        self.f_obs = masked_array(obs_cube[1].data * self.flux_unit, mask=np.isnan(obs_cube[1].data))
        self.f_err = masked_array(np.sqrt(obs_cube[2].data * self.flux_unit_err), mask=np.isnan(obs_cube[2].data))  # FIXME: CHECK!
        self.cube_shape = self.f_obs.shape

        # Model cube:
        model_cube = fits.open(sinopsis_directory + self.galaxy_id + '_modelcube.fits')[0]
        model_cube_nolines = fits.open(sinopsis_directory + self.galaxy_id + '_modelcube_nolines.fits')[0]

        self.f_syn = masked_array(model_cube.data * self.flux_unit, mask=model_cube.data == -999)
        self.f_syn_cont = masked_array(model_cube_nolines.data * self.flux_unit, mask=model_cube.data == -999)

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
            fname_fit_details = self.sinopsis_directory + fname_prefix + '.%0.4d_%0.4d.out' % (
            y + 1, x + 1)  # FIXME: 0- or 1-indexed?

            n_features = open(fname_fit_details, 'rt').read().splitlines()[0].split()[0]
            n_features = int(n_features)

            ssp_results = np.genfromtxt(fname_fit_details, skip_header=23)

            age, sfr = np.log10(ssp_results[:, 0]), ssp_results[:, 3]

            return age, sfr

            # fit_details_table = Table.read(fname_fit_details, format='ascii', header_start=1,
            #                                data_end=n_features+2)

            # ssp_results = Table.read(fname_fit_details, format='ascii', header_start=23,
            #                          data_start=24)

    def plot_map(self, sinopsis_property, show_plot=True, ax=None, custom_mask=None):

        sinplot.plot_sinopsis_map(self, sinopsis_property, ax=ax, custom_mask=custom_mask)

        if show_plot:
            plt.show()

    def plot_spectrum(self, x, y, plot_error=True, plot_legend=True, show_plot=True):

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        plt.figure()

        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y],
                         self.f_err[:, x, y], plot_error=plot_error, plot_legend=plot_legend, z=self.catalog['z'])

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

        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y],
                         self.f_err[:, x, y], ax=ax_spectrum, z=self.catalog['z'])

        sinplot.plot_residuals(self.wl, self.f_obs[:, x, y], self.f_syn_cont[:, x, y], ax=ax_residuals,
                               z=self.catalog['z'])

        # Plotting SFH:
        if self.config['sfh_type'] == 'ff':
            sinplot.plot_sfh(self.age_bin_center, self.sfh[:, x, y], ax=ax_sfh)

        if self.config['sfh_type'] == 'dexp':
            age, sfr = self.fit_details(x, y)

            ax_sfh.plot(age, sfr)

            ax_sfh.set_xlabel(r'$\log \, t \, \mathrm{[yr]}$', fontsize=12)
            ax_sfh.set_ylabel('SFR', fontsize=12)

        # FIXME: A problem with the units prevents reconstruction of the parametric SFHs
        #     t_initial = np.linspace(1e6, 1.4e9, 1000)
        #     t_late_burst = np.linspace(1e6, self.properties['Tb'][x, y], 1000)
        #     initial = utils.initial_burst(t_initial, t_u=1.4e9, n1=self.properties['n1'][x, y],
        #                                   tau_i=self.properties['tau_i'][x, y])
        #     late_burst = utils.late_burst(t_late_burst, m_b=self.properties['Mb'][x, y],
        #                                   t_b=self.properties['Tb'][x, y], n2=self.properties['n2'][x, y],
        #                                   tau_b=self.properties['tau_b'][x, y])
        #     ax_sfh.plot(np.log10(t_initial), initial, color='red')
        #     ax_sfh.plot(np.log10(t_late_burst), late_burst, color='blue')
        #     ax_sfh.set_ylabel('SFR', fontsize=12)
        #     ax_sfh.set_xlabel(r'$\log \, t \, \mathrm{[yr]}$', fontsize=12)
        #     ax_sfh.set_xlim(np.log10(1e6), np.log10(1.4e9))

        # A map to show the spaxel:
        # FIXME: Plotting map of total flux, did not think enough about this
        # FIXME: Inverting x and y !!! Have to find a better solution for this (numpy works in line, column not x, y)
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
            ax_residuals.annotate(annotations[i], xy=(0.82, 0.3 - i * 0.05), xycoords='figure fraction', fontsize=12)

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


class Sinopsis1D:
    """

    Reads SINOPSIS output for a catalog of 1d spectra into a python object.

    input:
    -----------
    sinopsis_directory: Directory with SINOPSIS files (must end with /)
        type: str

    """

    def __init__(self, sinopsis_directory):
        self.catalog_in = Table.read(sinopsis_directory + 'catalog.in', data_start=1, format='ascii')
        self.catalog_out = Table.read(sinopsis_directory + 'catalog.out', format='ascii')
        self.catalog_mag = Table.read(sinopsis_directory + 'catalog.mag', format='ascii')
        self.catalog_eqw = Table.read(sinopsis_directory + 'catalog.eqw', format='ascii')
        self.catalog_chi = Table.read(sinopsis_directory + 'catalog.chi', format='ascii')

        self.catalog_sfh = np.full(fill_value=np.nan, shape=(len(self.catalog_in), 12))

        for i in range(len(self.catalog_in)):
            self.catalog_sfh[i] = np.array([self.catalog_out['sfr_' + str(j)][i] for j in range(1, 13)])

        self.age_bins = np.genfromtxt(sinopsis_directory + 'catalog.log', skip_header=21, skip_footer=5)[:, 0]
        self.age_bins = np.append(0, self.age_bins)
        self.age_bins_4 = np.genfromtxt(sinopsis_directory + 'catalog.bin', skip_header=3)
        self.age_bins_4 = np.append(0, self.age_bins_4)


if __name__ == '__main__':
    import importlib

    importlib.reload(sinplot)
    importlib.reload(utils)

    sinopsis_cube = SinopsisCube('tests/test_run/')
    sinopsis_cube.plot_spectrum(54, 29, plot_error=False)
    sinopsis_cube.plot_map('SFR1')
    sinopsis_cube.plot_fit_complete(54, 29)

    # sinopsis_1d = Sinopsis1D('tests/test1d/')
