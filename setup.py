from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="pysinopsis",
    version="0.1.0",
    description="A set of python tools to work with SINOPSIS output.",
    # long_description=long_description,
    url="https://github.com/arielwrl/pysinopsis",
    author="Ariel Werle",
    author_email="ariel.werle@inaf.it",
    package_dir={"": "pysinopsis"},
    packages=find_packages(where="pysinopsis"),
    python_requires=">=3.2, <4",
    install_requires=["astropy>=5.2.2",
                      "matplotlib>=3.7.2",
                      "numpy>=1.17.4",
                      "scipy>=1.3.3",
                      "setuptools>=45.2.0",
                      "scikit-image>=0.0"], 
    include_package_data=True
)