"""

ariel@padova
01/06/2020

Tools for plotting SINOPSIS output.

"""

import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

# FIXME: Not sure of many of these labels
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
                   'lVwage': r'$\log \langle t_\star \rangle_V \mathrm{[yr]}$',
                   'lVwage_m': r'$\log \langle t_\star \rangle_V^\mathrm{min} \mathrm{[yr]}$',
                   'lVwage_M': r'$\log \langle t_\star \rangle_V^\mathrm{max} \mathrm{[yr]}$',
                   'lwage': r'$\log \langle t_\star \rangle_L \mathrm{[yr]}$',
                   'lwage_m': r'$\log \langle t_\star \rangle_L^\mathrm{min} \mathrm{[yr]}$',
                   'lwage_M': r'$\log \langle t_\star \rangle_L^\mathrm{max} \mathrm{[yr]}$',
                   'mwage': r'$\log \langle t_\star \rangle_M \mathrm{[yr]}$',
                   'mwage_m': r'$\log \langle t_\star \rangle_M^\mathrm{min} \mathrm{[yr]}$',
                   'mwage_M': r'$\log \langle t_\star \rangle_M^\mathrm{max} \mathrm{[yr]}$',
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
               'f_lambda': r'$F_\lambda \, \mathrm{[erg \, s^{-1} \, cm^{-2} \, \AA^{-1}]}$'}


def plot_fit(wl, f_obs, f_syn, f_syn_cont, f_err=None, ax=None, plot_error=True, plot_legend=True, flux_unit=1):
    """

    :param wl:
    :param f_obs:
    :param f_syn:
    :param f_syn_cont:
    :param f_err:
    :param ax:
    :param plot_error:
    :param plot_legend:
    :param flux_unit:

    """

    if plot_error is True and f_err is None:
        print('plot_error is set to True but no errors were passed')

    if ax is None:
        ax = plt.gca()

    ax.plot(wl, flux_unit * f_obs, color='k', lw=0.5, label='Observed Spectrum')
    ax.plot(wl, flux_unit * f_syn, color='blue', label='Synthetic Spectrum')
    ax.plot(wl, flux_unit * f_syn_cont, color='green', label='Synthetic Spectrum (Continuum)')

    if plot_error:
        ax.plot(wl, flux_unit * f_err, '--r', label='Error')

    if plot_legend:
        ax.legend()

    ax.set_xlabel(spec_labels['wl'])
    ax.set_ylabel(spec_labels['f_lambda'])


def plot_sinopsis_map(sinopsis_cube, sinopsis_property, cmap='magma_r', ax=None):

    if ax is None:
        ax = plt.gca()

    sinopsis_map = ax.imshow(sinopsis_cube.properties[sinopsis_property], cmap=cmap)

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    ax.figure.colorbar(mappable=sinopsis_map, label=sinopsis_labels[sinopsis_property])

