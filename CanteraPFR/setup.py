# -*- coding: utf-8 -*-
"""
Build Cython interface to CanteraPFR.

TODO
----
Automate system resolution (Nix vs Win).
"""

import os
from numpy import get_include
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

basepath = os.path.join(os.pardir, 'external', 'Nix')

sources = [
    'CanteraPFR.pyx',
    # os.path.join('src', 'AdiabaticPFR.cpp'),
    # os.path.join('src', 'HeatWallPFR.cpp'),
    # os.path.join('src', 'IsothermalPFR.cpp')
    ]

include_dirs = [
    # os.path.join(os.pardir, 'cantera', 'ext', 'sundials', 'include'),
    get_include(),
    os.path.join(basepath, 'include'),
    'include'
    ]

extra_objects = [
    os.path.join('lib', 'libCanteraPFR.a'),
    os.path.join(basepath, 'lib', 'libcantera.a')
    ]

extra_compile_args = []
extra_link_args = ['-lopenblas']
# '-lsundials_ida', '-lsundials_nvecserial',

setup(
    name = 'CanteraPFR',
    author = 'Walter Dal\'Maz Silva',
    author_email = 'waltermateriais@gmail.com',
    description = 'Plug-flow reactor models',
    url = 'https://github.com/waltermateriais/CanteraPFR/',
    license = 'UNLICENSE',
    version = '0.1.0',
    include_dirs = include_dirs,
    # FIXME add all packages here!
    # install_requires=['numpy>=1.11.1']
    ext_modules = cythonize(
        Extension(
            name = 'CanteraPFR',
            sources = sources,
            extra_compile_args = extra_compile_args,
            extra_objects = extra_objects,
            extra_link_args = extra_link_args,
            language = 'c++'
            ),
        build_dir = 'build',
        language_level = 3
    )
)
