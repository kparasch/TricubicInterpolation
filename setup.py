from setuptools import setup, Extension, find_packages
import sys
import numpy as np

if sys.version_info[0] < 3:
    extension_name = 'TricubicInterpolation/cTricubic/Tricubic2_c'
else:
    extension_name = 'TricubicInterpolation/cTricubic/Tricubic_c'

with open("README.md",  "r") as fh:
    long_description = fh.read()

module = Extension(extension_name,
                   ["TricubicInterpolation/cTricubic/source/coefs.c",
                    "TricubicInterpolation/cTricubic/source/coords.c",
                    "TricubicInterpolation/cTricubic/source/derivs.c",
                    "TricubicInterpolation/cTricubic/source/evaluate.c",
                    "TricubicInterpolation/cTricubic/source/pywrap.c",
                   ],
                   include_dirs=["TricubicInterpolation/cTricubic/include", np.get_include()],
                   extra_compile_args=["-std=c99"]
                  )


setup(
    name="TricubicInterpolation",
    version="1.1.0",
    author="Konstantinos Paraschou",
    author_email="konstantinos.paraschou@cern.ch",
    description="Tricubic Interpolation suited for Hamiltonian Particle Tracking written in C",
    long_description=long_description,
    url="https://github.com/kparasch/TricubicInterpolation",
    packages=find_packages(),
    keywords="Interpolation hamiltonian",
    ext_modules=[module]
)
