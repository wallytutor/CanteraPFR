# -*- coding: utf-8 -*-
"""
Build Cython interface to CanteraPFR.
"""

import os
import sys
from numpy import get_include
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


def get_platform():
    if 'linux' in sys.platform:
        platform_abbr = 'Nix'
    elif 'darwin' in sys.platform:
        platform_abbr = 'Nix'
    else:
        platform_abbr = 'Win'
    return platform_abbr


platform_abbr = get_platform()

# TODO use __file__ as initial path to avoid weird stuff!
basepath = os.path.abspath(os.path.join('external', platform_abbr))

include_dirs = [
    get_include(),
    os.path.join(basepath, 'include'),
    os.path.join('CanteraPFR', 'include')
    ]

extra_objects = [
    os.path.join('CanteraPFR', 'lib', 'libCanteraPFR.a'),
    os.path.join(basepath, 'lib', 'libcantera.a'),
    os.path.join(basepath, 'lib', 'libsundials_ida.a'),
    os.path.join(basepath, 'lib', 'libsundials_nvecserial.a')
    ]

extra_compile_args = []

if platform_abbr == 'Nix':
    extra_link_args = ['-lopenblas']
else:
    extra_link_args = []

# FIXME add all packages here!
# install_requires = ['numpy>=1.15.1']

ext_modules = []
ext_modules += [
    Extension(
        name = 'CanteraPFR.CanteraPFR',
        sources = [os.path.join('CanteraPFR', 'CanteraPFR.pyx')],
        extra_compile_args = extra_compile_args,
        extra_objects = extra_objects,
        extra_link_args = extra_link_args,
        language = 'c++'
        )
    ]
ext_modules += [
    Extension(
        name = 'CanteraPFR.CanteraAux',
        sources = [os.path.join('CanteraPFR', 'CanteraAux.pyx')],
        )
    ]

setup(
    name = 'CanteraPFR',
    packages = ['CanteraPFR'],
    author = 'Walter Dal\'Maz Silva',
    author_email = 'waltermateriais@gmail.com',
    description = 'Plug-flow reactor models',
    url = 'https://github.com/waltermateriais/CanteraPFR/',
    license = 'UNLICENSE',
    version = '0.1.1',
    include_dirs = include_dirs,
    ext_modules = cythonize(ext_modules,
                            build_dir = 'build',
                            language_level = 3
                            )
)
