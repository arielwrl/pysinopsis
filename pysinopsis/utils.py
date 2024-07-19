import numpy as np
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from skimage.measure import regionprops


def resampler(x_old, y_old, x_new):
    """
    Resamples y_old to match the new x_new grid.

    Parameters:
        x_old (array-like): Old x values.
        y_old (array-like): Old y values.
        x_new (array-like): New x values to interpolate onto.

    Returns:
        array-like: Resampled y values on the new grid.
    """

    interp = interp1d(x_old, y_old, bounds_error=False, fill_value=(0., 0.))

    y_new = interp(x_new)

    return y_new
    

def calc_sn(wl, f_obs, f_err, z, window_limits=(5500, 5700)):
    """
    Calculates the signal-to-noise ratio (SNR) within a specified wavelength window.

    Parameters:
        wl (array-like): Wavelength array.
        f_obs (array-like): Observed flux array.
        f_err (array-like): Flux error array.
        z (float): Redshift.
        window_limits (tuple): Lower and upper limits of the wavelength window (default: (5500, 5700)).

    Returns:
        float: Signal-to-noise ratio (SNR).
    """
 
    wl = wl / (1 + z)

    window = (wl > window_limits[0]) & (wl < window_limits[1])

    sn = np.mean(f_obs[window]) / np.mean(np.sqrt(f_err[window]))

    return sn


def calc_continuum_rms(wl, f_obs, f_syn, z, window_limits=(5500, 5700)):
    """
    Calculates the root mean square (RMS) deviation of the continuum flux between observed and synthetic spectra
    within a specified wavelength window.

    Parameters:
        wl (array-like): Wavelength array.
        f_obs (array-like): Observed flux array.
        f_syn (array-like): Synthetic flux array.
        z (float): Redshift.
        window_limits (tuple): Lower and upper limits of the wavelength window (default: (5500, 5700)).

    Returns:
        float: Root mean square (RMS) deviation of the continuum flux.
    """ 
    
    wl = wl / (1 + z)

    window = (wl > window_limits[0]) & (wl < window_limits[1])

    continuum_rms = np.mean(np.absolute(f_obs[window] - f_syn[window]))

    return continuum_rms


def get_uncertainty(sinopsis_cube, sinopsis_property, x, y):
    """
    Calculate the uncertainty in a given SINOPSIS property at a specific spaxel position.

    Parameters:
        sinopsis_cube (SinopsisCube): An instance of the SinopsisCube class.
        sinopsis_property (str): The SINOPSIS property for which uncertainty is to be calculated.
        x (int): x coordinate of the spaxel.
        y (int): y coordinate of the spaxel.

    Returns:
        tuple: A tuple containing the positive and negative uncertainties in the SINOPSIS property.
    """  

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
    Calculate the initial burst star formation rate.
    Equation (2) in page 29 of SINOPSIS manual, version 1.6.7

    Parameters:
        t (float): Time elapsed since the beginning of the burst.
        t_u (float): Duration of the burst.
        n1 (float): Power-law index.
        tau_i (float): Decay timescale of the burst.

    Returns:
        float: Initial burst star formation rate (SFR) at time t.
    """

    sfr = (((t_u - t) / t_u) ** n1) * np.exp(-((t_u - t) / (tau_i * t_u)))

    return sfr


def late_burst(t, m_b, t_b, n2, tau_b):
    """
    Calculate the late burst star formation rate (SFR) based on the SINOPSIS manual (version 1.6.7).

    Parameters:
        t (float): Time elapsed since the beginning of the burst.
        m_b (float): Mass of the burst (in Msun).
        t_b (float): Time at which the burst occurs.
        n2 (float): Power-law index.
        tau_b (float): Decay timescale of the burst.

    Returns:
        float: Late burst star formation rate (SFR) at time t.
    """

    sfr = m_b * (((t_b - t) / t_b) ** n2) * np.exp(-((t_b - t) / (tau_b * t_b)))

    return sfr


def gini(x):
    """
    Calculate the Gini coefficient for a given array (e.g the SFH).

    Parameters:
        x (array-like): Input array for which Gini coefficient will be calculated.

    Returns:
        float: Gini coefficient of the input array.
    """

    mad = np.abs(np.subtract.outer(x, x)).mean()
    rmad = mad/np.mean(x)
    g = 0.5 * rmad
    
    return g


def calc_mwage(age_bins_mid, sfrs, bin_widths):
    """
    Calculate the mass-weighted age.

    Parameters:
        age_bins_mid (array-like): Midpoints of age bins.
        sfrs (array-like): Star formation rates corresponding to each age bin.
        bin_widths (array-like): Widths of each age bin.

    Returns:
        float: Mass-weighted age of the stellar population.
    """

    mass_fractions = sfrs * bin_widths

    normalized_mass_fractions = mass_fractions / np.sum(mass_fractions)

    mwage = np.sum(age_bins_mid * normalized_mass_fractions)/np.sum(normalized_mass_fractions)

    return mwage


def cumulative_sfh(sfrs, bin_widths):
    """
    Calculate the normalized cumulative star formation history (cSFH).

    Parameters:
        sfrs (array-like): Star formation rates corresponding to each age bin.
        bin_widths (array-like): Widths of each age bin.

    Returns:
        array-like: Normalized cumulative SFH.
    """

    mass_fractions = sfrs * bin_widths

    normalized_mass_fractions = mass_fractions / np.sum(mass_fractions)

    cumulative_sfh = np.cumsum(normalized_mass_fractions[::-1])[::-1]

    return cumulative_sfh


def calc_t_x(age_bins, sfrs, bin_widths, target_fraction):
    """
    Calculate the age at which a certain fraction of the total stellar mass has formed.

    Parameters:
        age_bins (array-like): Ages corresponding to each bin.
        sfrs (array-like): Star formation rates corresponding to each age bin.
        bin_widths (array-like): Widths of each age bin.
        target_fraction (float): Target fraction of total stellar mass.

    Returns:
        float: Age at which the target fraction of total stellar mass has formed.
    """

    csfh = cumulative_sfh(sfrs, bin_widths)

    csfh = np.append(csfh, 0)

    interpolator = interp1d(age_bins, csfh)

    probe = np.linspace(age_bins[1], age_bins[-1], 40000)

    t_x = probe[np.argmin(np.absolute(interpolator(probe)-target_fraction))]

    return t_x


def calc_mass_formed(age_bins, sfrs, bin_widths, target_time):
    """
    Calculate the total stellar mass formed up to a certain age.

    Parameters:
        age_bins (array-like): Ages corresponding to each bin.
        sfrs (array-like): Star formation rates corresponding to each age bin.
        bin_widths (array-like): Widths of each age bin.
        target_time (float): Target age up to which the mass is formed.

    Returns:
        float: Total stellar mass formed up to the target age.
    """
    
    csfh = cumulative_sfh(sfrs, bin_widths)

    csfh = np.append(csfh, 0)

    interpolator = interp1d(age_bins, csfh)

    mass_formed = interpolator(target_time)

    return mass_formed


def calc_manual_ew(wl, flux, delta_wl, line):
    """
    Calculate the equivalent width (EW) of a spectral line manually.

    Parameters:
        wl (array-like): Wavelength array.
        flux (array-like): Flux array.
        delta_wl (float): Wavelength bin width.
        line (str): Name of the spectral line for which EW is calculated.
            Supported lines: 'Hd', 'Hb', 'Ha', 'Oii', 'Oiii'.

    Returns:
        float: Equivalent width of the specified spectral line.
    """

    if line == 'Hd':

        blue_cont = np.median(flux[(wl > 4076.4) & (wl < 4088.8)])
        red_cont = np.median(flux[(wl > 4117.2) & (wl < 4136.7)])

        blue_wl = np.mean(wl[(wl > 4076.4) & (wl < 4088.8)])
        red_wl = np.mean(wl[(wl > 4117.2) & (wl < 4136.7)])

        ew_range = (wl > 4090) & (wl < 4117)

    if line == 'Hb':

        blue_cont = np.median(flux[(wl > 4806.0) & (wl < 4826.0)])
        red_cont = np.median(flux[(wl > 4896.0) & (wl < 4918.0)])

        blue_wl = np.mean(wl[(wl > 4806.0) & (wl < 4826.0)])
        red_wl = np.mean(wl[(wl > 4896.0) & (wl < 4918.0)])

        ew_range = (wl > 4826.0) & (wl < 4896.0)

    if line == 'Ha':

        blue_cont = np.median(flux[(wl > 6505.0) & (wl < 6535.0)])
        red_cont = np.median(flux[(wl > 6595.0) & (wl < 6625.0)])

        blue_wl = np.mean(wl[(wl > 6505.0) & (wl < 6535.0)])
        red_wl = np.mean(wl[(wl > 6595.0) & (wl < 6625.0)])

        ew_range = (wl > 6553.0) & (wl < 6573.0)

    if line == 'Oii':

        blue_cont = np.median(flux[(wl > 3670.0) & (wl < 3718.0)])
        red_cont = np.median(flux[(wl > 3738.0) & (wl < 3774.0)])

        blue_wl = np.mean(wl[(wl > 3670.0) & (wl < 3717.0)])
        red_wl = np.mean(wl[(wl > 3738.0) & (wl < 3774.0)])

        ew_range = (wl > 3717.0) & (wl < 3737.0)
    
    if line == 'Oiii':

        blue_cont = np.median(flux[(wl > 4985.5) & (wl < 4998)])
        red_cont = np.median(flux[(wl > 5018.0) & (wl < 5036.0)])

        blue_wl = np.mean(wl[(wl > 4985.5) & (wl < 4998)])
        red_wl = np.mean(wl[(wl > 5018.0) & (wl < 5036.0)])

        ew_range = (wl > 4997.0) & (wl < 5017.0)
    
    cont_slope = (red_cont - blue_cont) / (red_wl - blue_wl)
    cont_intercept = blue_cont - blue_wl * cont_slope

    cont_level = wl * cont_slope + cont_intercept

    integrated_hdelta = np.trapz(1 - (flux[ew_range] / cont_level[ew_range]), dx=delta_wl)

    return integrated_hdelta


def box_filter(wl_in, flux_in, box_width=16):
    """
    Apply a box filter to smooth the input spectrum.

    Parameters:
        wl_in (array-like): Input wavelength array.
        flux_in (array-like): Input flux array.
        box_width (int): Width of the box filter window.

    Returns:
        tuple: Tuple containing the smoothed wavelength array and smoothed flux array.
    """

    wl_out = np.array([np.mean(wl_in[i:i+box_width]) for i in range(len(wl_in)-box_width)])
    flux_out = np.array([np.mean(flux_in[i:i + box_width]) for i in range(len(flux_in) - box_width)])

    return wl_out, flux_out


def luminosity_distance(z, h0=70, omega0=0.3, in_cm=False):
    """
    Calculate the luminosity distance using the specified cosmological parameters.

    Parameters:
        z (float): Redshift.
        h0 (float, optional): Hubble constant in km/s/Mpc. Default is 70.
        omega0 (float, optional): Omega0 parameter. Default is 0.3.
        in_cm (bool, optional): If True, return luminosity distance in cm. Default is False.

    Returns:
        float: Luminosity distance, default in Mpc.
    """

    cosmo = FlatLambdaCDM(h0, omega0)   

    dl = cosmo.luminosity_distance(z)

    if in_cm:
        dl = dl.to(u.cm)

    return dl.value


def calc_center_of_mass(image, label_image):
    """
    Calculate the center of mass of an object in an image based on its labeled image.

    Parameters:
        image (ndarray): Input image.
        label_image (ndarray): Labeled image containing the object.

    Returns:
        tuple: Coordinates of the center of mass (row, column).
    """

    # label_image = label_image.astype(int)
    image_properties = regionprops(label_image, intensity_image=image)[0]
    center_of_mass = image_properties.centroid_weighted

    return center_of_mass

