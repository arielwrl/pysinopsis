"""

ariel@padova
15/06/2020

Script to plot SINOPSIS fits of each individual spaxel

"""

import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from pysinopsis.output import SinopsisCube
import os

sns.set_style('whitegrid')

galaxy_id = 'A370_01'
sinopsis_dir = 'C:/Users/ariel/Workspace/GASP/High-z/SINOPSIS/A370_01/'
plot_dir = '../tests/test_plots'
plot_format = 'png'

sinopsis_cube = SinopsisCube(sinopsis_dir)

if galaxy_id not in os.listdir(plot_dir):
    print(plot_dir+galaxy_id, 'not found, creating directory...')
    os.makedirs(plot_dir+galaxy_id)

for i, j in itertools.product(range(sinopsis_cube.cube_shape[1]), range(sinopsis_cube.cube_shape[2])):

    print(i, j)

    if sinopsis_cube.invalid_spaxel(i, j):
        print('Skipping invalid spaxel')
        continue
    else:
        plot_file_name = plot_dir + galaxy_id + '/' + '%s_%i_%i' % (galaxy_id, i, j) + '.' + plot_format
        print('Plotting spaxel to', plot_file_name)

    sinopsis_cube.plot_fit_complete(i, j, out_file_name=plot_file_name, out_format=plot_format, show_plot=False)
