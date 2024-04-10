import numpy as np
import matplotlib.pyplot as plt
import os
import re


class ModelSet:
    """
    
    Creates a ModelSet object with the following attributes:

    models_dir: Directory containing the model files.
    file_list: List of files in the models_dir.
    model_list: List of model files with the '.dat' extension.
    metallicities: Array of metallicities parsed from the model filenames.
    log_ages: Array of logarithmic ages parsed from the model filenames.
    ages: Array of ages calculated from the logarithmic ages.
    wl: Wavelength array extracted from one of the model files.
    model_spectra: Array of model spectra within the specified wavelength range.

    """

    def __init__(self, models_dir='./', wl_range=[1000, 10000]):
        """
        Initializes a ModelSet object.

        Parameters:
            models_dir (str): Directory containing the model files.
            wl_range (list): Wavelength range of interest.
            """

        self.models_dir = models_dir

        self.file_list = os.listdir(models_dir)
        self.model_list = [file for file in self.file_list if file.split('.')[-1] == 'dat']
       
        self.metallicities = np.array([float(re.split('Z|_|t|.dat', model_name)[3]) for model_name in self.model_list])
        self.log_ages = np.array([float(re.split('Z|_|t|.dat', model_name)[5]) for model_name in self.model_list])
        self.ages = np.array([10**log_age for log_age in self.log_ages])

        self.wl = np.genfromtxt(self.models_dir + self.model_list[1]).transpose()[0]

        wl_flag = (self.wl > wl_range[0]) & (self.wl < wl_range[1])

        self.wl = self.wl[wl_flag]
        self.model_spectra = np.array([np.genfromtxt(self.models_dir + model_name).transpose()[1][wl_flag]
                                       for model_name in self.model_list])


if __name__ == '__main__':

    cb2020_dir = 'C:/Users/ariel/Workspace/sinopsis/data/ssp/ssp_cb2020/'

    sinopsis_models = ModelSet(cb2020_dir)

    for i in [0, 400, 1000]:
        plt.plot(sinopsis_models.wl, sinopsis_models.model_spectra[i], 
                 label='$Z=$' + str(sinopsis_models.metallicities[i]) + ', $\log\,t=$' +
                 str(sinopsis_models.log_ages[i]))
    
    plt.xlim(2000, 7000)
    plt.legend()