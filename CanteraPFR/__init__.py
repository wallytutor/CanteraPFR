# -*- coding: utf-8 -*-
from .ct_pfr import *
from .ct_aux import *
from .ct_test import *


def set_env():
    import os
    import cantera
    absdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    absdir = os.path.join(absdir, 'data')
    os.environ['CANTERA_DATA'] = f'{absdir}:$CANTERA_DATA'
    cantera.add_directory(absdir)
    cantera.suppress_thermo_warnings()


def test_all():
    set_env()
    test_dimensionless()
    test_filter_mechanism()
    test_analyse_graph()
    test_pypfr()


set_env()


__version__ = '0.1.2'
