import numpy as np
import matplotlib.pyplot as plt
from pysinopsis.models import ModelSet
from pysinopsis.utils import luminosity_distance
from pysinopsis.synphot import synflux, synmag
from astropy.table import Table


cb2020_dir = 'C:/Users/ariel/Workspace/sinopsis/data/ssp/ssp_cb2020/'
f275w_filter = 'C:/Users/ariel/Workspace/pysinopsis/pysinopsis/data/filters/HST_WFC3_UVIS2.F275W.dat'

sinopsis_models = ModelSet(cb2020_dir, wl_range=[1500, 3500])

redshift = 0.056
distance_cm = luminosity_distance(z=redshift, in_cm=True)

wl = sinopsis_models.wl * (1 + redshift)

model_fluxes = 1e30 * 10 ** sinopsis_models.model_spectra / (4 * np.pi * distance_cm ** 2)

model_fluxes_f275w = np.array([synflux(sinopsis_models.wl, model_fluxes[i], f275w_filter) for i
                                   in range(len(model_fluxes))])
model_magnitudes_f275w = np.array([synmag(sinopsis_models.wl, model_fluxes[i], f275w_filter) for i
                                   in range(len(model_fluxes))])

models_table = Table()
models_table['age'] = sinopsis_models.ages / 1e9
models_table['metallicity'] = sinopsis_models.metallicities / 0.017
models_table['f275w_flux'] = model_fluxes_f275w
models_table['f275w_mag'] = model_magnitudes_f275w

models_table = models_table[models_table['age'] < 1.1]

models_table.write('C:/Users/ariel/Workspace/pysinopsis/pysinopsis/data/models_table.rst',
                   format='ascii.rst', overwrite=True)