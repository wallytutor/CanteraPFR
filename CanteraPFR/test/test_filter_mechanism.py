# -*- coding: utf-8 -*-

from CanteraPFR import filter_mechanism

gas = filter_mechanism('gri30.cti', ['C2H2', 'H2', 'C2H4', 'CH4'])
print(gas.report())

gas = filter_mechanism('gri30.cti', ['C2H2', 'H2', 'C2H4', 'CH4'], write=True,
                       idname='gas', output='test_filter_mechanism.xml',
                       overwrite=True)
print(gas.report())
