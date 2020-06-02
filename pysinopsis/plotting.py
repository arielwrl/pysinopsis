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
                   'Av_y_m': r'$$',
                   'Av_y_M': r'$$',
                   'Av': r'$A_V$',
                   'Av_m': r'$$',
                   'Av_M': r'$$',
                   'SFR1': r'$SFR_1$',
                   'SFR1_m': r'$$',
                   'SFR1_M': r'$$',
                   'SFR2': r'$SFR_2$',
                   'SFR2_m': r'$$',
                   'SFR2_M': r'$$',
                   'SFR3': r'$SFR_3$',
                   'SFR3_m': r'$$',
                   'SFR3_M': r'$$',
                   'SFR4': r'$SFR_4$',
                   'SFR4_m': r'$$',
                   'SFR4_M': r'$$',
                   'Mb1_3': r'$$',
                   'Mb2_3': r'$$',
                   'Mb3_3': r'$$',
                   'Mb4_3': r'$$',
                   'Mb1_2': r'$$',
                   'Mb1_m_2': r'$$',
                   'Mb1_M_2': r'$$',
                   'Mb2_2': r'$$',
                   'Mb2_m_2': r'$$',
                   'Mb2_M_2': r'$$',
                   'Mb3_2': r'$$',
                   'Mb3_m_2': r'$$',
                   'Mb3_M_2': r'$$',
                   'Mb4_2': r'$$',
                   'Mb4_m_2': r'$$',
                   'Mb4_M_2': r'$$',
                   'Mb1_1': r'$$',
                   'Mb2_1': r'$$',
                   'Mb3_1': r'$$',
                   'Mb4_1': r'$$',
                   'AMass3': r'$$',
                   'AMass3_m': r'$$',
                   'AMass3_M': r'$$',
                   'TotMass3': r'$$',
                   'AMass2': r'$$',
                   'AMass2_m': r'$$',
                   'AMass2_M': r'$$',
                   'TotMass2': r'$$',
                   'AMass1': r'$$',
                   'AMass1_m': r'$$',
                   'AMass1_M': r'$$',
                   'TotMass1': r'$$',
                   'lVwage': r'$$',
                   'lVwage_m': r'$$',
                   'lVwage_M': r'$$',
                   'lwage': r'$$',
                   'lwage_m': r'$$',
                   'lwage_M': r'$$',
                   'mwage': r'$$',
                   'mwage_m': r'$$',
                   'mwage_M': r'$$',
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

    map = ax.imshow(sinopsis_cube.properties[sinopsis_property], cmap=cmap)

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    ax.figure.colorbar(mappable=map, label=sinopsis_labels[sinopsis_property])

