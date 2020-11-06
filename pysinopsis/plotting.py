"""

ariel@padova
01/06/2020

Tools for plotting SINOPSIS output.

"""

import numpy as np
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
                   'lVwage': r'$\langle t_\star \rangle_V \mathrm{[Gyr]}$',
                   'lVwage_m': r'$\langle t_\star \rangle_V^\mathrm{min} \mathrm{[Gyr]}$',
                   'lVwage_M': r'$\langle t_\star \rangle_V^\mathrm{max} \mathrm{[Gyr]}$',
                   'lwage': r'$\langle t_\star \rangle_L \mathrm{[Gyr]}$',
                   'lwage_m': r'$\langle t_\star \rangle_L^\mathrm{min} \mathrm{[Gyr]}$',
                   'lwage_M': r'$\langle t_\star \rangle_L^\mathrm{max} \mathrm{[Gyr]}$',
                   'mwage': r'$\langle t_\star \rangle_M \mathrm{[Gyr]}$',
                   'mwage_m': r'$\langle t_\star \rangle_M^\mathrm{min} \mathrm{[Gyr]}$',
                   'mwage_M': r'$\langle t_\star \rangle_M^\mathrm{max} \mathrm{[Gyr]}$',
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

    FIXME: Doc me

    :param z:
    :param syn_lw:
    :param obs_lw:
    :param syn_cont_color:
    :param syn_color:
    :param obs_color:
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

    if ax is None:
        ax = plt.gca()

    ax.plot(wl / (1+z), (f_obs-f_syn_cont) * (1+z) / f_syn_cont, color=res_color, lw=res_lw)

    ax.set_ylabel(spec_labels['res'], fontsize=12)
    ax.set_xlabel(spec_labels['wl'], fontsize=12)

    ax.set_ylim(-0.35, 0.35)

    ax.hlines(y=0, xmin=wl[0] / (1+z), xmax=wl[-1] / (1+z), lw=2, zorder=15, linestyles='dashed')


def plot_sinopsis_map(sinopsis_cube, sinopsis_property, cmap='magma_r', ax=None, custom_mask=None):

    if ax is None:
        ax = plt.gca()

    if custom_mask is not None:
        sinopsis_map = ax.imshow(np.ma.masked_array(sinopsis_cube.properties[sinopsis_property], mask=custom_mask),
                                 cmap=cmap, origin='lower')
    else:
        sinopsis_map = ax.imshow(sinopsis_cube.properties[sinopsis_property], cmap=cmap, origin='lower')

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    cb = plt.colorbar(mappable=sinopsis_map, ax=ax)
    cb.set_label(sinopsis_labels[sinopsis_property], fontsize=20)


def plot_sfh(age_bin_center, sfh_array, ax=None):

    if ax is None:
        ax = plt.gca()

    ax.scatter(np.log10(age_bin_center), sfh_array)

    ax.set_ylabel('SFR', fontsize=12)
    ax.set_xlabel(r'$\log \, t \, \mathrm{[yr]}$', fontsize=12)

