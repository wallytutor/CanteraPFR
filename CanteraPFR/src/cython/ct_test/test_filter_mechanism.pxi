# -*- coding: utf-8 -*-


def test_filter_mechanism():
    from .ct_aux import filter_mechanism

    species = ['C2H2', 'H2', 'C2H4', 'CH4']
    gas = filter_mechanism('gri30.cti', species)
    print(gas.report())

    gas = filter_mechanism('gri30.cti', species, write=True, idname='gas',
                           output='test_filter_mechanism.xml', overwrite=True)
    print(gas.report())
