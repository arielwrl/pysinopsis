import matplotlib.pyplot as plt
from pysinopsis.models import ModelSet

cb2020_dir = 'C:/Users/ariel/Workspace/sinopsis/data/ssp/ssp_cb2020/'

sinopsis_models = ModelSet(cb2020_dir)

for i in [0, 400, 1000]:
    plt.plot(sinopsis_models.wl, sinopsis_models.model_spectra[i], 
                label='$Z=$' + str(sinopsis_models.metallicities[i]) + ', $\log\,t=$' +
                str(sinopsis_models.log_ages[i]))

plt.xlim(2000, 7000)
plt.legend()