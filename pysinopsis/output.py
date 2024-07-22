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
    Reads configuration settings from a SINOPSIS configuration file.

    Parameters:
        sinopsis_dir (str, optional): The directory containing the SINOPSIS configuration file. 
            Defaults to './'.

    Returns:
        dict: A dictionary containing the extracted configuration settings.
            - 'input_catalog' (str): The input catalog file name.
            - 'input_type' (str): The type of input data (e.g., 'cube', '1D Spectrum').
            - 'sinopsis_dir' (str): The directory where SINOPSIS was run.
            - 'sfh_type' (str): The type of star formation history (e.g., 'ff', 'lognormal').

    Raises:
        FileNotFoundError: If the SINOPSIS configuration file is not found.
    """
 
    config_file = np.array(open(sinopsis_dir + 'config.sin', 'rb').readlines()).astype(str)

    config_dict = {'input_catalog': config_file[12].split()[-1],
                   'input_type': config_file[19].split()[-1],
                   'sinopsis_dir': sinopsis_dir,
                   'sfh_type': config_file[62].split()[-1],
                   }

    if config_dict['input_type'] == 'cube':
        config_dict['galaxy_id'] = config_dict['input_catalog'].split('.')[0]

    return config_dict


def read_sinopsis_catalog(sinopsis_dir, input_catalog_file, input_type):
    """
    Read SINOPSIS catalog information.

    Parameters:
        sinopsis_dir (str): Directory containing SINOPSIS data.
        input_catalog_file (str): Filename of the input catalog file.
        input_type (str): Type of input data ('cube', 'spec', etc.).

    Returns:
        dict: Dictionary containing SINOPSIS catalog information.
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


def read_results_cube(output_cube_file, cube_type):
    """
    Read results data cube into a dictionary of masked arrays.

    Parameters:
        output_cube_file (str): Filename of the output cube file.
        cube_type (str): Type of cube ('dexp', 'GASP', etc.).

    Returns:
        tuple: A tuple containing header information and properties.
            The header information is a dictionary containing selected header keys.
            The properties is an OrderedDict containing properties extracted from the cube.
    """

    output_cube = fits.open(output_cube_file)[0]

    header_info_keys = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'EXTEND', 'CRVAL1', 'CRVAL2',
                        'CD1_1', 'CD2_1', 'CD1_2', 'CD2_2', 'CRPIX1', 'CRPIX2', 'CTYPE1', 'CTYPE2', 'OBJECT']

    header_info = OrderedDict()
    for key in header_info_keys:
        header_info[key] = output_cube.header[key]

    properties = OrderedDict()

    if cube_type == 'dexp':
        len_cube = 74
    if cube_type == 'GASP':
        len_cube = 77
    else:
        len_cube = 89

    for i in range(len_cube):
        plane_key = 'PLANE%0.2d' % i

        properties[output_cube.header[plane_key]] = masked_array(output_cube.data[i], mask=output_cube.data[i] == -999)

    return header_info, properties


def read_eqw_cube(eqw_cube_file):
    """
    Read SINOPSIS equivalent width (EQW) cube from a FITS file.

    Parameters:
        eqw_cube_file (str): Filename of the EQW cube FITS file.

    Returns:
        OrderedDict: A dictionary containing EQW data for each plane.
            The keys are EQW plane names, and the values are masked arrays representing EQW data.
    """

    eqw_cube = fits.open(eqw_cube_file)[0]

    eqws = OrderedDict()
    for i in range(eqw_cube.data.shape[0]):
        plane_key = 'PLANE%0.2d' % i  # FIXME: Hard-coded!

        eqws[eqw_cube.header[plane_key]] = masked_array(eqw_cube.data[i], mask=eqw_cube.data[i] == -999)

    return eqws


def read_mag_cube(mag_cube_file):
    """
    Read SINOPSIS magnitudes cube from a FITS file.

    Parameters:
        mag_cube_file (str): Filename of the magnitudes cube FITS file.

    Returns:
        OrderedDict: A dictionary containing magnitudes.
            The keys are magnitudes plane names, and the values are masked arrays representing magnitudes data.
    """

    mag_cube = fits.open(mag_cube_file)[0]

    mags = OrderedDict()
    for i in range(mag_cube.data.shape[0]):
        plane_key = 'PLANE%0.2d' % i  # FIXME: Hard-coded!

        mags[mag_cube.header[plane_key]] = masked_array(mag_cube.data[i], mask=mag_cube.data[i] == -999)

    return mags


def read_sfh_file(fname):
    """
    Read star formation history (SFH) data from a file.

    Parameters:
        fname (str): The filename of the SFH file.

    Returns:
        tuple: A tuple containing three numpy arrays:
            - age (numpy.ndarray): Array of ages.
            - sfr (numpy.ndarray): Array of star formation rates (SFRs).
            - ebv (numpy.ndarray): Array of extinction values (E(B-V)).
   
    Notes:
        The SFH file is expected to have a specific format with a header row and data rows.
        The function automatically detects the header row and skips it when reading the data.
    """

    file = open(fname, 'rb').readlines()
    for i in range(len(file)):
        if np.any(np.array(file[i].split())==b'Age'):
            skip_lines=i+1

    ssp_results = np.genfromtxt(fname, skip_header=skip_lines)

    age, sfr, ebv = ssp_results[:, 0], ssp_results[:, 3], ssp_results[:, 2]

    return age, sfr, ebv


def read_fit_results(fname, skip_header=1, skip_footer=12):
    """
    Read fitting results from a file.

    Parameters:
        fname (str): The filename of the file containing fitting results.
        skip_header (int): Number of header lines to skip when reading the file. Default is 1.
        skip_footer (int): Number of footer lines to skip when reading the file. Default is 12.

    Returns:
        dict: A dictionary containing fitting results as numpy arrays with keys corresponding 
        to parameter names.
    """

    results = np.genfromtxt(fname, skip_header=skip_header, skip_footer=skip_footer)

    results_dict = {'Lambc' : results[:, 0],
                    'Deltl' : results[:, 1],
                    'Sigma' : results[:, 2],
                    'Weight' : results[:, 3],
                    'Obs' : results[:, 5],
                    'Model' : results[:, 5],
                    'Chi' : results[:, 6]}

    return results_dict


def read_obs_cube(fname):

        obs_cube = fits.open(fname)

        obs_header = {'primary': obs_cube[0].header,
                      'flux': obs_cube[1].header,
                      'error': obs_cube[2].header}

        try:
            flux_unit = 10 ** float(obs_header['flux']['BUNIT'].split()[0].split('**')[1].split('(')[-1].split(')')[0])
        except Exception:
            print('Cannot read flux unit from ', obs_header['flux']['BUNIT'], ', trying to split()[0] instead:',
                  obs_header['flux']['BUNIT'].split()[0])
            flux_unit = float(obs_header['flux']['BUNIT'].split()[0])
        try:
            flux_unit_err = 10 ** float(obs_header['error']['BUNIT'].split()[0].split('**')[1].split('(')[-1].split(')')[0])
        except Exception:
            print('Cannot read flux unit from ', obs_header['error']['BUNIT'], ', trying to split()[0] instead:',
                  obs_header['error']['BUNIT'].split()[0])
            flux_unit_err = float(obs_header['error']['BUNIT'].split()[0])

        wl = obs_cube[1].header['CRVAL3'] + obs_cube[1].header['CD3_3'] * np.arange(obs_cube[1].header['NAXIS3'])
        f_obs = masked_array(obs_cube[1].data * flux_unit, mask=np.isnan(obs_cube[1].data))
        f_err = masked_array(np.sqrt(obs_cube[2].data * flux_unit_err), mask=np.isnan(obs_cube[2].data))  # FIXME: CHECK!
        cube_shape = f_obs.shape

        return {'obs_header': obs_header, 'flux_unit': flux_unit, 'flux_unit_err': flux_unit_err,
                'wl': wl, 'f_obs': f_obs, 'f_err': f_err, 'cube_shape': cube_shape}


class SinopsisCube:
    """
    
    Represents the full SINOPSIS output cubes and provides methods to manipulate and visualize the data.

    """

    def __init__(self, sinopsis_directory='./', memory_saving=False, verbose=False):
        """
        Initializes a SinopsisCube object.

        For very large input data cubes (e.g 4Gb), `memory_saving=True` is recommended.

        Parameters:
            sinopsis_directory (str): Directory containing SINOPSIS data.
            memory_saving (bool): Flag indicating whether to load certain data into memory or not. Default is False.
            verbose (bool): Flag indicating whether to print verbose information during initialization. Default is False.
        """

        config = read_config(sinopsis_directory)
        catalog = read_sinopsis_catalog(sinopsis_directory, config['input_catalog'], input_type='cube')

        self.config = config
        self.catalog = catalog

        self.galaxy_id = config['galaxy_id']
        self.obs_file = catalog['obs_file']
        self.sinopsis_directory = sinopsis_directory

        # SINOPSIS results:
        self.header_info, self.properties = read_results_cube(sinopsis_directory + self.galaxy_id + '_out.fits',
                                                              cube_type=self.config['sfh_type'])
        
        if verbose:
            print('Results correctly read')

        # Try to read velocity dispersion
        try:
            self.velocity_dispersion = fits.open(sinopsis_directory + self.galaxy_id + '_sig_abs_mask.fits')[0].data
        except Exception:
            self.velocity_dispersion = np.full_like(self.properties['mwage'], 0)

        # Ages:
        # TODO: Add proper support for other SFH types
        self.age_bins = np.genfromtxt(sinopsis_directory + self.galaxy_id + '.log', skip_header=22,
                                        skip_footer=6)[:, 0]  # FIXME: Double-check me
        self.n_ages = len(self.age_bins)
        self.age_bins = np.append(0, self.age_bins)
        self.age_bin_width = np.array([self.age_bins[i+1] - self.age_bins[i] for i in range(self.n_ages)])
        self.age_bins_4 = np.genfromtxt(sinopsis_directory + self.config['input_catalog'].split('.')[0] + '.bin', skip_header=3)
        self.age_bins_4 = np.append(0, self.age_bins_4)

        self.age_bin_center = np.array([(self.age_bins[i] + self.age_bins[i + 1]) / 2 for i in range(self.n_ages)])

        self.sfh = np.array([self.properties['sfr_' + str(i)] for i in range(1, self.n_ages+1)])
        self.sfh = np.ma.masked_array(self.sfh, mask=self.sfh == -999)

        self.integrated_sfh = self.sfh.sum(axis=(1, 2))

        # Equivalend widths:
        if ~memory_saving:
            self.eqw = read_eqw_cube(sinopsis_directory + self.galaxy_id + '_eqw.fits')

        # Magnitudes:
        if ~memory_saving:
            self.mag = read_mag_cube(sinopsis_directory + self.galaxy_id + '_mag.fits')

        # Mask:
        self.mask = fits.open(sinopsis_directory + self.galaxy_id + '_fitmask.fits')[0].data

        # Observed cube:
        obs_cube = read_obs_cube(self.sinopsis_directory + self.obs_file)
        if verbose:
            print('Finished reading observed cube')

        self.obs_header = obs_cube['obs_header']
        self.flux_unit = obs_cube['flux_unit']
        self.flux_unit_err = obs_cube['flux_unit_err']
        self.wl = obs_cube['wl']
        self.f_obs = obs_cube['f_obs']
        self.f_err = obs_cube['f_err']
        self.cube_shape = obs_cube['cube_shape']

        # Model cube:
        model_cube = fits.open(sinopsis_directory + self.galaxy_id + '_modelcube.fits')[0]
        model_cube_nolines = fits.open(sinopsis_directory + self.galaxy_id + '_modelcube_nolines.fits')[0]

        if verbose:
            print('Model cube read')

        self.f_syn = masked_array(model_cube.data * self.flux_unit, mask=model_cube.data == -999)
        self.f_syn_cont = masked_array(model_cube_nolines.data * self.flux_unit, mask=model_cube.data == -999)

        if verbose:
            print('f_syn Calculated')

        # Emission only cube:
        if ~memory_saving:
            self.emission_only = self.f_obs - self.f_syn_cont

        if ~memory_saving:
            self.emission_only_model = self.f_syn - self.f_syn_cont 

        if verbose:
            print('Emission-only calculated')


    def invalid_spaxel(self, x, y):
        """
        Checks if the given spaxel is invalid.

        Parameters:
            x (int): x-coordinate of the spaxel.
            y (int): y-coordinate of the spaxel.

        Returns:
            bool: True if the spaxel is invalid, False otherwise.
        """

        if np.any(self.f_syn[:, x, y].mask):
            return True
        else:
            return False


    def fit_details(self, x, y):
        """
        Retrieves the fitting details for a specific spaxel.

        Parameters:
            x (int): x-coordinate of the spaxel.
            y (int): y-coordinate of the spaxel.

        Returns:
            tuple: A tuple containing arrays of ages and SFR (star formation rate) values.
        """

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        else:

            fname_prefix = self.obs_file.split('.')[0]
            fname_fit_details = self.sinopsis_directory + fname_prefix + '.%0.4d_%0.4d.out' % (
            y + 1, x + 1)  

            age, sfr = read_sfh_file(fname_fit_details)

            return age, sfr


    def plot_map(self, sinopsis_property, show_plot=True, ax=None, custom_mask=None, cmap='magma_r'):
        """
        Plots a map of a SINOPSIS property.

        Parameters:
            sinopsis_property (str): Property whose map is to be plotted.
            show_plot (bool): Flag indicating whether to display the plot. Default is True.
            ax (matplotlib.axes.Axes): Axes object to plot on. If None, creates a new plot. Default is None.
            custom_mask (numpy.ndarray): Custom mask for masking certain areas in the map. Default is None.
            cmap (str): Name of the colormap to use. Default is 'magma_r'.
        """

        sinplot.plot_sinopsis_map(self, sinopsis_property, ax=ax, custom_mask=custom_mask, cmap=cmap)

        if show_plot:
            plt.show()
    

    def get_center_of_mass(self, sinopsis_property, custom_mask=None):
        """
        Calculates the center of mass for a given SINOPSIS property.

        Parameters:
            sinopsis_property (str): Property for which the center of mass is to be calculated.
            custom_mask (numpy.ndarray): Custom mask for masking certain areas. Default is None.

        Returns:
            tuple: A tuple containing the x and y coordinates of the center of mass.
        """

        if custom_mask is not None:
            label_image = custom_mask
        else:
            label_image = self.mask == 1
            
        center_of_mass = utils.calc_center_of_mass(self.properties[sinopsis_property], label_image)
            
        return center_of_mass


    def plot_spectrum(self, x, y, plot_error=True, plot_legend=True, show_plot=True, ax=None, obs_color='k',
                      syn_color='g'):
        """
        Plots the spectrum of a spaxel.

        Parameters:
            x (int): x-coordinate of the spaxel.
            y (int): y-coordinate of the spaxel.
            plot_error (bool): Flag indicating whether to plot error spectrum. Default is True.
            plot_legend (bool): Flag indicating whether to include a legend in the plot. Default is True.
            show_plot (bool): Flag indicating whether to display the plot. Default is True.
            ax (matplotlib.axes.Axes): Axes object to plot on. If None, creates a new plot. Default is None.
            obs_color (str): Color for observed spectrum. Default is 'k'.
            syn_color (str): Color for synthetic spectrum. Default is 'g'.
        """

        if ax is None:
            plt.figure()
            ax = plt.gca()

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], self.f_syn[:, x, y], self.f_syn_cont[:, x, y],
                         self.f_err[:, x, y], plot_error=plot_error, plot_legend=plot_legend, z=self.catalog['z'],
                         ax=ax, syn_color=syn_color, obs_color=obs_color)

        if show_plot:
            plt.show()


    def plot_fit_complete(self, x, y, figsize=(13.25, 8.5), out_file_name=None, out_format='png',
                          show_plot=True):
        """
        Plots a complete diagnostics of a fit for a given spaxel.

        Parameters:
            x (int): x-coordinate of the spaxel.
            y (int): y-coordinate of the spaxel.
            figsize (tuple): Size of the figure (width, height) in inches. Default is (13.25, 8.5).
            out_file_name (str): Output file name if saving the plot. If None, doesn't save. Default is None.
            out_format (str): Format for the output file. Default is 'png'.
            show_plot (bool): Flag indicating whether to display the plot. Default is True.
        """

        if self.invalid_spaxel(x, y):
            print('>>> Masked spaxel!')

        fig = plt.figure(figsize=figsize)

        gs = gridspec.GridSpec(5, 5)

        ax_spectrum = plt.subplot(gs[0:2, 0:5])
        ax_residuals = plt.subplot(gs[2, 0:5], sharex=ax_spectrum)
        ax_sfh = plt.subplot(gs[3:5, 0:2])
        ax_map = plt.subplot(gs[3:5, 2:4])

        # Plotting fit and residuals:
        fsyn = self.f_syn[:, x, y]
        fsyn_cont = self.f_syn_cont[:, x, y]

        sinplot.plot_fit(self.wl, self.f_obs[:, x, y], fsyn, fsyn_cont,
                         self.f_err[:, x, y], ax=ax_spectrum, z=self.catalog['z'])

        sinplot.plot_residuals(self.wl, self.f_obs[:, x, y], fsyn_cont, ax=ax_residuals,
                               z=self.catalog['z'])

        # Plotting SFH:
        if self.config['sfh_type'] == 'ff':
            sinplot.plot_sfh(self.age_bin_center, self.sfh[:, x, y], ax=ax_sfh)

        if self.config['sfh_type'] == 'dexp':
            age, sfr = self.fit_details(x, y)

            ax_sfh.plot(age, sfr)

            ax_sfh.set_xlabel(r'$\log \, t \, \mathrm{[yr]}$', fontsize=12)
            ax_sfh.set_ylabel('SFR', fontsize=12)

        # A map to show the spaxel:
        reference_map = ax_map.imshow(np.sum(self.f_obs, axis=0), cmap='Blues', origin='lower')
        fig.colorbar(mappable=reference_map, ax=ax_map, label=r'$\log \, \Sigma \, F_\lambda$')
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

    Provides minimal support for 1D fits.

    """

    def __init__(self, sinopsis_directory):
        """
        Initializes a Sinopsis1D object.

        Parameters:
            sinopsis_directory (str): The directory containing SINOPSIS files.
        """
        
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

