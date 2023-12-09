"""

ariel@prague
04/12/2023

Reads and organizes stellar population models in the standard SINOPSIS format

"""


import numpy as np
import matplotlib.pyplot as plt
import os
import re


class ModelSet:

    def __init__(self, models_dir='./'):

        self.models_dir = models_dir

        self.file_list = os.listdir(models_dir)
        self.model_list = [file for file in self.file_list if file.split('.')[-1] == 'dat']
        
        self.metallicities = [float(re.split('Z|_|t|.dat', model_name)[3]) for model_name in self.model_list]
        self.log_ages = [float(re.split('Z|_|t|.dat', model_name)[5]) for model_name in self.model_list]
        self.ages = [10**log_age for log_age in self.log_ages]

        self.wl = np.genfromtxt(self.models_dir + self.model_list[1]).transpose()[0]

        self.model_spectra = [np.genfromtxt(self.models_dir + model_name).transpose()[1] for model_name in self.model_list]


if __name__ == '__main__':

    cb2020_dir = 'C:/Users/ariel/Workspace/sinopsis/data/ssp/ssp_cb2020/'

    sinopsis_models = ModelSet(cb2020_dir)

    for i in [0, 400, 1000]:
        plt.plot(sinopsis_models.wl, sinopsis_models.model_spectra[i], 
                 label='$Z=$'+str(sinopsis_models.metallicities[i])+', $\log\,t=$'+
                 str(sinopsis_models.log_ages[i]))
    
    plt.xlim(2000, 7000)
    plt.legend()