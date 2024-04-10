# pysinopsis
###### A set of python tools to work with SINOPSIS output.

[SINOPSIS](https://www.irya.unam.mx/gente/j.fritz/JFhp/SINOPSIS.html) (SImulatiNg OPtical Spectra wIth Stellar population models) is a spectrophotometric fitting code which aims at reproducing the combined spectral and broad-band photometric data of galaxies.

The code has been extensively used to analyze galaxy data cubes, and it generates a large amount of data products that can be difficult to organize.

Here we provide a pythonic solution that automates most of the tasks dealing with the analysis and validation of SINOPSIS data products. The package is focused in the use case of data cubes, but minimal support for single spectra is also provided by the `Sinopsis1D()` class. 

## Instalation

To install `pysinopsis` you can run `pip install pysinopsis` or clone the repository and install from `setup.py`.

## Quickstart

The core class of the package is `SinopsisCube()`, and it can be initialized as follows:

```
from pysinopsis.output import SinopsisCube

sinopsis_dir = 'path/to/sinopsis/directory/'
sinopsis_cube = SinopsisCube(sinopsis_dir)
```

Where `sinopsis_dir` is the directory where sinopsis was run and contains all the input, output and config information for the data cube of a given galaxy. By default it is set to the current directory.

A more complete guide can be found in [this notebook](https://github.com/arielwrl/pysinopsis/blob/master/examples/quickstart.ipynb).

