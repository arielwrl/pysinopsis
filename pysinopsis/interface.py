import os
from pysinopsis.output import read_config, read_sinopsis_catalog


def smart_remove_file(file_name):
    """
    Remove a file if it exists, otherwise print a message indicating that the file was not found.

    Parameters:
        file_name (str): Name of the file to be removed.

    """
    
    try:
        os.remove(file_name)
    except FileNotFoundError:
        print('Did not find', file_name)


def clear_sinopsis_run(sinopsis_directory='./'):
    """
    Clears files generated in a SINOPSIS runs for datacubes

    Parameters:
        sinopsis_directory (str): Directory path where SINOPSIS files are located. Default is the current directory.
   
    """

    config = read_config(sinopsis_directory)

    galaxy_id = config['galaxy_id']

    smart_remove_file(sinopsis_directory + galaxy_id + '.log')
    smart_remove_file(sinopsis_directory + galaxy_id + '.bin')
    smart_remove_file(sinopsis_directory + galaxy_id + '_eqw.fits')
    smart_remove_file(sinopsis_directory + galaxy_id + '_fitmask.fits')
    smart_remove_file(sinopsis_directory + galaxy_id + '_fluxcont.fits')
    smart_remove_file(sinopsis_directory + galaxy_id + '_mag.fits')
    smart_remove_file(sinopsis_directory + galaxy_id + '_modelcube.fits')
    smart_remove_file(sinopsis_directory + galaxy_id + '_modelcube_nolines.fits')
    smart_remove_file(sinopsis_directory + galaxy_id + '_out.fits')

    print('>>> Removed all files!')


def write_cat_file(cat_fname, cube_fname, mask_fnames, redshift):
    """
    Write catalog file containing cube filename, mask filenames, and redshift.

    Parameters:
        cat_fname (str): Catalog file name.
        cube_fname (str): Filename of the data cube.
        mask_fnames (list): List of mask filenames.
        redshift (float): Redshift value.

    """

    cat_file = open(cat_fname, 'wb')

    cat_file.write(str.encode(cube_fname+'\n'))

    for i in range(len(mask_fnames)):
        cat_file.write(str.encode(mask_fnames[i] + ' '))
    cat_file.write(b'\n')

    cat_file.write(str.encode(str(redshift)))

    cat_file.close()






