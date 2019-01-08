# -*- coding: utf-8 -*-


def test_filter_mechanism():
    from .ct_aux import filter_mechanism
    
    spc = ['C2H2', 'H2', 'C2H4', 'CH4']
    gas = filter_mechanism('gri30.cti', spc)
    print(gas.report())

    gas = filter_mechanism('gri30.cti', spc, write=True, idname='gas',
                           output='test_filter_mechanism.xml', overwrite=True)
    print(gas.report())
