# -*- coding: utf-8 -*-
""" Setup script for building CanteraPFR extensions. """

from pathlib import Path
from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np

ROOT_DIR = Path(__file__).resolve().parent

ext_modules = [
    Extension(
        name='cantera_pfr.ct_pfr',
        sources=[
            'src/cantera_pfr/ct_pfr.pyx',
            'src/extension/AdiabaticPFR.cpp',
            'src/extension/HeatWallPFR.cpp',
            'src/extension/IsothermalPFR.cpp',
        ],
        include_dirs=[
            np.get_include(),
            str(ROOT_DIR / 'src' / 'extension'),
            '/usr/local/include',
            '/usr/include/eigen3',
            '/usr/lib/x86_64-linux-gnu/openmpi/include',
        ],
        libraries=[
            'cantera_shared',
            'sundials_ida',
            'sundials_nvecserial',
            'openblas',
            'lapack',
            'pthread',
        ],
        library_dirs=[
            '/usr/local/lib',
            '/usr/lib/x86_64-linux-gnu/openmpi/lib',
        ],
        extra_compile_args=['-std=c++17', '-O3'],
        extra_link_args=['-Wl,-rpath,/usr/local/lib'],
        language='c++',
    )
]

setup(
    ext_modules=cythonize(
        ext_modules,
        language_level=3,
    )
)
