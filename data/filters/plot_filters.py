import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from toolbox import plot_tools
import bagpipes as pipes

sns.set_style('ticks')


def plot_filters(paths_to_filters, filter_names, colors = sns.diverging_palette(220, 20, n=5)):

    fig = plt.figure(figsize=(5.5, 4.25))
    
    ages = [0.005]

    for i in range(1):

        exp = {}
        exp["age"] = ages[i]
        # exp["tau"] = 10
        exp['metallicity'] = 1
        exp['massformed'] = 10

        dust = {}
        dust["type"] = "Cardelli"
        dust["Av"] = 0.1

        nebular = {}
        nebular["logU"] = -2.75
        nebular["velshift"] = 0.0

        dust["eta"] = 1

        model_components = {}
        model_components["redshift"] = 0.044631
        model_components["burst"] = exp
        model_components["dust"] = dust
        model_components["veldisp"] = 75
        model_components["nebular"] = nebular

        # Wavelength array to plot spectra.
        wl = np.arange(2000, 9900, 2)

        # Creating a model galaxy
        model = pipes.model_galaxy(model_components, filt_list=paths, phot_units='ergscma', spec_wavs=wl)

        plt.plot(wl, model.spectrum[:, 1] / np.max(model.spectrum[:, 1]),
                c='k', lw=0.3, zorder=0)


    for i in range(len(paths_to_filters)):
        filter = paths_to_filters[i]


        if i in [0, 2, 3]:

            wl, transmission = np.genfromtxt(filter).transpose()
            plt.plot(wl, transmission, color=colors[i], lw=1, label=filter_names[i])
            
        else:
            
            wl, transmission = np.genfromtxt(filter).transpose()
            plt.plot(wl, transmission, color=colors[i], lw=1, label=filter_names[i], ls='dashed')
            
        
            # plt.fill_between(wl, np.zeros_like(transmission), transmission,
            #                 color=colors[i], alpha=0.1)


    plt.legend(frameon=False, fontsize=11, loc=0)

    plt.ylim(0, 0.3)

    plt.xlabel(r'$\lambda \mathrm{[\AA]}$', fontsize=16)
    plt.ylabel(r'$T_\lambda$', fontsize=16)

    plt.tick_params(axis='both', labelsize=11)

    sns.despine()

    fig.subplots_adjust(top=0.965,
                        bottom=0.136,
                        left=0.125,
                        right=0.968,
                        hspace=0.2,
                        wspace=0.2)


filters = ['F275W', 'F336W', 'F606W', 'F680N', 'F814W']
paths = ['c:/Users/ariel/Workspace/GASP/HST/Data/filters/HST_WFC3_UVIS2.' + filter_name + '.dat' for filter_name in filters]

colors = ['mediumvioletred', '#0058b0', 'indigo', 'goldenrod', '#990147']


plot_filters(paths, filters, colors=colors)
# plt.savefig('filters.jpg', dpi=300)
# plt.savefig('filters.pdf')

plt.xlim(2000, 9900)

plt.show()