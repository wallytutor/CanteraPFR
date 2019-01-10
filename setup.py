# -*- coding: utf-8 -*-
"""
Build Cython interface to CanteraPFR.
"""

import os
import sys
from numpy import get_include
from setuptools import setup
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
absdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
basepath = os.path.join(absdir, 'external', platform_abbr)

if not os.path.exists(basepath):
    print('Assuming default path Cantera and Sundials')
    basepath = ''

install_requires = [
    'numpy>=1.11.0',
    'pandas>=0.23.4',
    'matplotlib>=2.2.0',
    'networkx>=2.1',
    'cantera>=2.4.0'
]

include_dirs = [
    get_include(),
    os.path.join(basepath, 'include'),
    os.path.join(absdir, 'CanteraPFR', 'include')
    ]

extra_objects = [
    os.path.join(absdir, 'CanteraPFR', 'lib', 'libCanteraPFR.a'),
    os.path.join(basepath, 'lib', 'libcantera.a'),
    os.path.join(basepath, 'lib', 'libsundials_ida.a'),
    os.path.join(basepath, 'lib', 'libsundials_nvecserial.a')
    ]

extra_compile_args = []

if platform_abbr == 'Nix':
    extra_link_args = ['-lopenblas']
else:
    extra_link_args = []

ext_modules = []

ext_modules += [
    Extension(
        name = 'CanteraPFR.ct_pfr',
        sources = [os.path.join('CanteraPFR', 'ct_pfr.pyx')],
        extra_compile_args = extra_compile_args,
        extra_objects = extra_objects,
        extra_link_args = extra_link_args,
        language = 'c++'
        )
    ]

ext_modules += [
    Extension(
        name = 'CanteraPFR.ct_aux',
        sources = [os.path.join('CanteraPFR', 'ct_aux.pyx')],
        )
    ]

ext_modules += [
    Extension(
        name = 'CanteraPFR.ct_graph',
        sources = [os.path.join('CanteraPFR', 'ct_graph.pyx')],
        )
    ]

ext_modules += [
    Extension(
        name = 'CanteraPFR.ct_test',
        sources = [os.path.join('CanteraPFR', 'ct_test.pyx')],
        )
    ]

ext_modules = cythonize(ext_modules, build_dir='build', language_level=3)

setup(
    name = 'CanteraPFR',
    packages = ['CanteraPFR'],
    include_package_data = True,
    install_requires = install_requires,
    include_dirs = include_dirs,
    ext_modules = ext_modules,
    author = 'Walter Dal\'Maz Silva',
    author_email = 'waltermateriais@gmail.com',
    description = 'Plug-flow reactor models',
    keywords = 'kinetics reactor chemistry cantera transport',
    url = 'https://github.com/waltermateriais/CanteraPFR/',
    license = 'UNLICENSE',
    version = '0.1.2',
    zip_safe = False
)
