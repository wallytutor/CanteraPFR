# -*- coding: utf-8 -*-
from .ct_pfr import *
from .ct_aux import *
from .ct_test import *

def set_env():
    import os
    absdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    os.environ['CANTERA_DATA'] = os.path.join(absdir, 'data')
    # FIXME this is not holding as expected!
    # ct.add_directory(os.path.join(absdir, 'data'))
    # ct.suppress_thermo_warnings()


def test_all():
    test_dimensionless()
    test_filter_mechanism()
    test_pypfr()


set_env()


__version__ = '0.1.2'
