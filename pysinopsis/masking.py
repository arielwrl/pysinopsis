import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.patches import Circle
from astropy.io import fits

plt.ion()  # God forgive me for hardcoding this, such an astronomy thing to do


def read_image(image_file):
    """
    Read image data from a FITS file.

    Args:
        image_file (str): Filename of the FITS file containing the image data.

    Returns:
        numpy.ndarray: Image data.
    """
 

    image = fits.open(image_file)[1].data

    return image


def plot_image(image, vmin=0, vmax=15):
    """
    Plot an image.

    Parameters:
        image (numpy.ndarray): The image data to be plotted.
        vmin (float, optional): Minimum value for the color scale. Default is 0.
        vmax (float, optional): Maximum value for the color scale. Default is 15.

    """

    plt.imshow(image.transpose(), vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.show()


def draw_circle(xc, yc, radius, axis=None):
    """
    Draw a circle on a matplotlib axis.

    Parameters:
        xc (float): x-coordinate of the center of the circle.
        yc (float): y-coordinate of the center of the circle.
        radius (float): Radius of the circle.
        axis (matplotlib.axis.Axis, optional): The axis on which to draw the circle. If None, the current axis is used. Default is None.

    Returns:
        matplotlib.patches.Circle: The drawn circle object.
    """

    if axis is None:
        axis = plt.gca()

    circle = Circle((xc, yc), radius=radius)

    axis.add_patch(circle)

    return circle


def create_mask(xc, yc, radius, image):
    """
    Create a circular mask on an image centered at (xc, yc) with a given radius.

    Parameters:
        xc (int): x-coordinate of the center of the circle.
        yc (int): y-coordinate of the center of the circle.
        radius (int): Radius of the circle.
        image (numpy.ndarray): Input image.

    Returns:
        numpy.ndarray: The spatial mask.
    """

    n_x = image.shape[0]
    n_y = image.shape[1]

    x, y = np.ogrid[-yc:n_x - yc, -xc:n_y - xc]

    masked = x ** 2 + y ** 2 > radius ** 2

    spatial_mask = np.zeros((n_x, n_y))
    spatial_mask[masked] = 1.

    plt.gca().imshow(np.ma.masked_array(spatial_mask, mask=spatial_mask == 1).transpose(), zorder=10, alpha=0.5)

    return spatial_mask


def write_mask_fits(spatial_mask, out_fname):
    """
    Write a mask array to a FITS file.

    Parameters:
        spatial_mask (numpy.ndarray): The spatial mask array.
        out_fname (str): Output FITS file name.
    """

    hdu = fits.PrimaryHDU(spatial_mask)
    hdu_list = fits.HDUList([hdu])

    hdu_list.writeto(out_fname, overwrite=True)
