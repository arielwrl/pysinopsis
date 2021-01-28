"""

ariel@bonporti

Tools for creating custom masks for datacubes.

Make sure matplotlib is in interactive mode.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.patches import Circle
from astropy.io import fits

plt.ion()  # God fogive me for hardcoding this, such an astronomy thing to do


def read_image(image_file):

    image = fits.open(image_file)[1].data

    return image


def plot_image(image, vmin=0, vmax=15):

    plt.imshow(image.transpose(), vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.show()


def draw_circle(xc, yc, radius, axis=None):

    if axis is None:
        axis = plt.gca()

    circle = Circle((xc, yc), radius=radius)

    axis.add_patch(circle)

    return circle


def create_mask(xc, yc, radius, image):

    n_x = image.shape[0]
    n_y = image.shape[1]

    x, y = np.ogrid[-yc:n_x - yc, -xc:n_y - xc]

    masked = x ** 2 + y ** 2 > radius ** 2

    spatial_mask = np.zeros((n_x, n_y))
    spatial_mask[masked] = 1.

    plt.gca().imshow(np.ma.masked_array(spatial_mask, mask=spatial_mask==1).transpose(), zorder=10)

    return spatial_mask


def write_mask_fits(spatial_mask, out_fname):

    hdu = fits.PrimaryHDU(spatial_mask)
    hdu_list = fits.HDUList([hdu])

    hdu_list.writeto(out_fname, overwrite=True)


if __name__=='__main__':

    image_fname='/home/ariel/Workspace/pysinopsis/tests/masking_example/MACS1206_11_DATACUBE_FINAL_v1_SDSS_r.fits'

    image = read_image(image_fname)

    plot_image(image, vmin=0, vmax=10)

    # mask_patch = draw_circle(28, 23.5, 4)

    mask = create_mask(23.5, 28, 5, image)

    write_mask_fits(mask, '/home/ariel/Workspace/pysinopsis/tests/masking_example/MACS1206_11_mask.fits')
