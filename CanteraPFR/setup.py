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

if 'linux' in sys.platform:
    platform_abbr = 'Nix'
elif 'darwin' in sys.platform:
    platform_abbr = 'Nix'
else:
    platform_abbr = 'Win'

# TODO use __file__ as initial path to avoid weird stuff!
basepath = os.path.abspath(os.path.join(os.pardir, 'external', platform_abbr))

sources = ['CanteraPFR.pyx']
include_dirs = [get_include(), os.path.join(basepath, 'include'), 'include']

extra_objects = [
    os.path.join('lib', 'libCanteraPFR.a'),
    os.path.join(basepath, 'lib', 'libcantera.a'),
    os.path.join(basepath, 'lib', 'libsundials_ida.a'),
    os.path.join(basepath, 'lib', 'libsundials_nvecserial.a')
    ]

extra_compile_args = []
extra_link_args = ['-lopenblas']

# FIXME add all packages here!
# install_requires = ['numpy>=1.15.1']

setup(
    name = 'CanteraPFR',
    author = 'Walter Dal\'Maz Silva',
    author_email = 'waltermateriais@gmail.com',
    description = 'Plug-flow reactor models',
    url = 'https://github.com/waltermateriais/CanteraPFR/',
    license = 'UNLICENSE',
    version = '0.1.0',
    include_dirs = include_dirs,
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
