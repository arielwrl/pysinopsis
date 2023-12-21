import numpy as np
from scipy.interpolate import interp1d


def resampler(x_old, y_old, x_new):

    interp = interp1d(x_old, y_old, bounds_error = False
    , fill_value = (0.,0.))

    y_new = interp(x_new)
    
    return y_new
    

def synflux(wl, flux, filter_curve):
    
    if type(filter_curve) is str: 
        wl_filter, T = np.genfromtxt(filter_curve).transpose()
    else:
        wl_filter, T = filter_curve[0], filter_curve[1]
    
    wl_new = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5)
    
    T = resampler(wl_filter, T, wl_new)
    flux = resampler(wl, flux, wl_new)
    
    synflux = np.trapz(flux * T * wl_new, dx=1)
    synflux /= np.trapz(wl_new * T, dx=1)

    return synflux


def synmag(wl, flux, filter_curve, error=None, flag=None, badpix_tolerance=0.25, 
           interpolate_bad_pixels=False):

    if type(filter_curve) is str: 
        wl_filter, T = np.genfromtxt(filter_curve).transpose()
    else:
        wl_filter, T = filter_curve[0], filter_curve[1]

    if flag is not None:
        filter_range = (wl > wl_filter[0]) & (wl < wl_filter[-1])
        if flag[filter_range].sum() > badpix_tolerance * len(flag[filter_range]):
            badpix = True
        else:
            badpix = False

    wl_new = np.arange(np.round(wl_filter[0]) - 5, np.round(wl_filter[-1]) + 5)

    if interpolate_bad_pixels:
        T = resampler(wl_filter, T, wl_new)
        flux = resampler(wl[~flag], flux[~flag], wl_new)
        if error is not None:
            error = resampler(wl[~flag], error[~flag], wl_new)
    else: 
        T = resampler(wl_filter, T, wl_new)
        flux = resampler(wl, flux, wl_new)
        if error is not None:
            error = resampler(wl, error, wl_new)

    m_ab = -2.5 * np.log10(np.trapz(flux * T * wl_new, dx=1) / np.trapz(T / wl_new, dx=1)) - 2.41

    if error is not None:
        m_ab_error = 1.0857 * np.sqrt(np.sum(T**2 * error**2 * wl_new ** 2)) / np.sum(flux * T * wl_new)
        if flag is not None:
            return m_ab, m_ab_error, badpix
        else:
            return m_ab, m_ab_error
    else:
        return m_ab


def pivot_wavelength(filter_curve):

    if type(filter_curve) is str: 
        wl_filter, T = np.genfromtxt(filter_curve).transpose()
    else:
        wl_filter, T = filter_curve[0], filter_curve[1]
    
    pivot_wl = np.trapz(T * wl_filter, dx=1) / np.trapz(T * (wl_filter**-1), dx=1)
    pivot_wl = np.sqrt(pivot_wl)
    
    return pivot_wl


def effective_wavelength(wl, spectrum, filter_curve):

    if type(filter_curve) is str: 
        wl_filter, T = np.genfromtxt(filter_file).transpose()
    else:
        wl_filter, T = filter_curve[0], filter_curve[1]
    
    wl_filter_new = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5,1)
    
    T_filter_new = resampler(wl_filter, T_filter, wl_filter_new)
    spectrum_new = resampler(wl, spectrum, wl_filter_new)
    
    effective_wl = np.trapz(wl_filter_new**2 * spectrum_new * T_filter_new, dx=1)
    effective_wl /= np.trapz(spectrum_new * T_filter_new * wl_filter_new, dx=1)
    
    return effective_wl
    
    
    
    
    
    
