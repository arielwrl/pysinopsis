import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True) 

sinopsis_labels = {'Dl': r'$D_L \mathrm{[Mpc]}$',
                   'z': r'$z$',
                   'redchi': r'$\chi^2/N$',
                   'Z': r'$Z/Z_\odot$',
                   'nrun': r'$n_{\mathrm{Run}}$',
                   'neqw': r'$n_{EW}$',
                   'Av_y': r'$A_V^Y$',
                   'Av_y_m': r'$A_V^Y \mathrm{min}$',
                   'Av_y_M': r'$A_V^Y \mathrm{max}$',
                   'Av': r'$A_V$',
                   'Av_m': r'$A_V \mathrm{min}$',
                   'Av_M': r'$A_V \mathrm{max}$',
                   'SFR1': r'$SFR_1$',
                   'SFR1_m': r'$SFR_1^\mathrm{min}$',
                   'SFR1_M': r'$SFR_1^\mathrm{max}$',
                   'SFR2': r'$SFR_2$',
                   'SFR2_m': r'$SFR_2^\mathrm{min}$',
                   'SFR2_M': r'$SFR_2^\mathrm{max}$',
                   'SFR3': r'$SFR_3$',
                   'SFR3_m': r'$SFR_3^\mathrm{min}$',
                   'SFR3_M': r'$SFR_3^\mathrm{max}$',
                   'SFR4': r'$SFR_4$',
                   'SFR4_m': r'$SFR_4^\mathrm{min}$',
                   'SFR4_M': r'$SFR_4^\mathrm{max}$',
                   'Mb1_3': r'$\mu_1$',
                   'Mb2_3': r'$\mu_2$',
                   'Mb3_3': r'$\mu_3$',
                   'Mb4_3': r'$\mu_4$',
                   'Mb1_2': r'$\mu_1$',
                   'Mb1_m_2': r'$\mu_1^\mathrm{min}$',
                   'Mb1_M_2': r'$\mu_1^\mathrm{max}$',
                   'Mb2_2': r'$\mu_2$',
                   'Mb2_m_2': r'$\mu_2^\mathrm{min}$',
                   'Mb2_M_2': r'$\mu_2^\mathrm{max}$',
                   'Mb3_2': r'$\mu_3$',
                   'Mb3_m_2': r'$\mu_3^\mathrm{min}$',
                   'Mb3_M_2': r'$\mu_3^\mathrm{max}$',
                   'Mb4_2': r'$\mu_4$',
                   'Mb4_m_2': r'$\mu_4^\mathrm{min}$',
                   'Mb4_M_2': r'$\mu_4^\mathrm{max}$',
                   'Mb1_1': r'$\mu_1$',
                   'Mb2_1': r'$\mu_2$',
                   'Mb3_1': r'$\mu_3$',
                   'Mb4_1': r'$\mu_4 $',
                   'AMass3': r'$M_\star/M_\odot$',
                   'AMass3_m': r'$M_\star/M_\odot ^\mathrm{min}$',
                   'AMass3_M': r'$M_\star/M_\odot ^\mathrm{max}$',
                   'TotMass3': r'$M_\star/M_\odot$',
                   'AMass2': r'$M_\star/M_\odot$',
                   'AMass2_m': r'$M_\star/M_\odot ^\mathrm{min}$',
                   'AMass2_M': r'$M_\star/M_\odot ^\mathrm{max}$',
                   'TotMass2': r'$M_\star/M_\odot$',
                   'AMass1': r'$M_\star/M_\odot$',
                   'AMass1_m': r'$M_\star/M_\odot ^\mathrm{min}$',
                   'AMass1_M': r'$M_\star/M_\odot ^\mathrm{max}$',
                   'TotMass1': r'$M_\star/M_\odot$',
                   'lVwage': r'$\langle \log \, t_\star \rangle_V \mathrm{[yr]}$',
                   'lVwage_m': r'$\langle \log \, t_\star \rangle_V^\mathrm{min} \mathrm{[yr]}$',
                   'lVwage_M': r'$\langle \log \, t_\star \rangle_V^\mathrm{max} \mathrm{[yr]}$',
                   'lwage': r'$\langle \log \, t_\star \rangle_L \mathrm{[yr]}$',
                   'lwage_m': r'$\langle \log \, t_\star \rangle_L^\mathrm{min} \mathrm{[yr]}$',
                   'lwage_M': r'$\langle \log \, t_\star \rangle_L^\mathrm{max} \mathrm{[yr]}$',
                   'mwage': r'$\langle \log \, t_\star \rangle_M \mathrm{[yr]}$',
                   'mwage_m': r'$\langle \log \, t_\star \rangle_M^\mathrm{min} \mathrm{[yr]}$',
                   'mwage_M': r'$\langle \log \, t_\star \rangle_M^\mathrm{max} \mathrm{[yr]}$',
                   'sfr_1': r'$SFR_1$',
                   'sfr_2': r'$SFR_2$',
                   'sfr_3': r'$SFR_3$',
                   'sfr_4': r'$SFR_4$',
                   'sfr_5': r'$SFR_5$',
                   'sfr_6': r'$SFR_6$',
                   'sfr_7': r'$SFR_7$',
                   'sfr_8': r'$SFR_8$',
                   'sfr_9': r'$SFR_9$',
                   'sfr_10': r'$SFR_10$',
                   'sfr_11': r'$SFR_11$'}

spec_labels = {'wl': r'$\lambda \, \mathrm{[\AA]}$',
               'f_lambda': r'$F_\lambda \, \mathrm{[erg \, s^{-1} \, cm^{-2} \, \AA^{-1}]}$',
               'res': r'$\left(O_\lambda - M_\lambda^C \right)/M_\lambda^C$'}


def plot_fit(wl, f_obs, f_syn, f_syn_cont, f_err=None, ax=None, plot_error=True, plot_legend=True, flux_unit=1,
             flux_unit_err=1, obs_color='k', syn_color='b', syn_cont_color='g', obs_lw=0.5, syn_lw=1.5, z=0):
    """
    Plot observed and synthetic spectra.

    Parameters:
      wl (array-like): Wavelength array.
      f_obs (array-like): Observed flux array.
      f_syn (array-like): Synthetic flux array.
      f_syn_cont (array-like): Synthetic continuum flux array.
      f_err (array-like, optional): Flux error array. Default is None.
      ax (matplotlib.axes.Axes, optional): Axes object for plotting. If None, creates new. Default is None.
      plot_error (bool, optional): Whether to plot error spectrum. Default is True.
      plot_legend (bool, optional): Whether to plot the legend. Default is True.
      flux_unit (float, optional): Unit conversion factor for flux. Default is 1.
      flux_unit_err (float, optional): Unit conversion factor for flux error. Default is 1.
      obs_color (str, optional): Color for observed spectrum. Default is 'k' (black).
      syn_color (str, optional): Color for synthetic spectrum. Default is 'b' (blue).
      syn_cont_color (str, optional): Color for synthetic continuum spectrum. Default is 'g' (green).
      obs_lw (float, optional): Linewidth for observed spectrum. Default is 0.5.
      syn_lw (float, optional): Linewidth for synthetic spectra. Default is 1.5.
      z (float, optional): Redshift. Default is 0.
    """

    if plot_error is True and f_err is None:
        print('plot_error is set to True but no errors were passed')

    if ax is None:
        ax = plt.gca()

    ax.plot(wl / (1+z), flux_unit * f_obs * (1+z), color=obs_color, lw=obs_lw, label='Observed Spectrum')
    ax.plot(wl / (1+z), flux_unit * f_syn * (1+z), color=syn_color, lw=syn_lw, label='Synthetic Spectrum')
    ax.plot(wl / (1+z), flux_unit * f_syn_cont * (1+z), color=syn_cont_color, lw=syn_lw,
            label='Synthetic Spectrum (Continuum)')

    if plot_error:
        ax.plot(wl / (1+z), flux_unit_err * f_err * (1+z), '--r', label='Error')

    if plot_legend:
        ax.legend()

    ax.set_xlabel(spec_labels['wl'], fontsize=12)
    ax.set_ylabel(spec_labels['f_lambda'], fontsize=12)


def plot_residuals(wl, f_obs, f_syn_cont, res_color='g', res_lw=0.5, z=0, ax=None):
    """
    Plot residuals between observed and synthetic continuum spectra.

    Parameters:
      wl (array-like): Wavelength array.
      f_obs (array-like): Observed flux array.
      f_syn_cont (array-like): Synthetic continuum flux array.
      res_color (str, optional): Color for residuals plot. Default is 'g' (green).
      res_lw (float, optional): Linewidth for residuals plot. Default is 0.5.
      z (float, optional): Redshift. Default is 0.
      ax (matplotlib.axes.Axes, optional): Axes object for plotting. If None, uses the current Axes. Default is None.
    """

    if ax is None:
        ax = plt.gca()

    ax.plot(wl / (1+z), (f_obs-f_syn_cont) * (1+z) / f_syn_cont, color=res_color, lw=res_lw)

    ax.set_ylabel(spec_labels['res'], fontsize=12)
    ax.set_xlabel(spec_labels['wl'], fontsize=12)

    ax.hlines(y=0, xmin=wl[0] / (1+z), xmax=wl[-1] / (1+z), lw=2, zorder=15, linestyles='dashed')


def plot_sinopsis_map(sinopsis_cube, sinopsis_property, cmap='magma_r', ax=None, custom_mask=None):
    """
    Plot a map of a SINOPSIS property.

    Parameters:
      sinopsis_cube (SinopsisCube): The SINOPSIS cube containing the property data.
      sinopsis_property (str): The property to plot.
      cmap (str, optional): The colormap for the plot. Default is 'magma_r'.
      ax (matplotlib.axes.Axes, optional): Axes object for plotting. If None, uses the current Axes. Default is None.
      custom_mask (array-like, optional): Custom mask to apply to the property data. Default is None.
    """
    
    if ax is None:
        ax = plt.gca()

    if custom_mask is not None:
        sinopsis_map = ax.imshow(np.ma.masked_array(sinopsis_cube.properties[sinopsis_property], mask=custom_mask),
                                 cmap=cmap, origin='lower')
    else:
        sinopsis_map = ax.imshow(sinopsis_cube.properties[sinopsis_property], cmap=cmap, origin='lower')

    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('y', fontsize=16)

    cb = plt.colorbar(mappable=sinopsis_map, ax=ax)
    cb.set_label(sinopsis_labels[sinopsis_property], fontsize=20)


def plot_sfh(age_bin_center, sfh_array, ax=None, sfh_color='steelblue'):
    """
    Plot the star formation history (SFH).

    Parameters:
      age_bin_center (array-like): Array containing the central ages of the bins.
      sfh_array (array-like): Array containing the star formation rates (SFRs).
      ax (matplotlib.axes.Axes, optional): Axes object for plotting. If None, uses the current Axes. Default is None.
    """
    
    if ax is None:
        ax = plt.gca()

    ax.scatter(np.log10(age_bin_center), sfh_array, color=sfh_color)
    ax.plot(np.log10(age_bin_center), sfh_array, color=sfh_color)

    ax.set_ylabel('SFR', fontsize=12)
    ax.set_xlabel(r'$\log \, t \, \mathrm{[yr]}$', fontsize=12)

